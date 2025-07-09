classdef PropofolSystem

    % an automated intravenous propofol infusion control system
    % using a three-compartment mamillary pkpd model. 
    % implements an extended kalman filter and barrier function
    % controllers. 

    properties

        % plant settings
        A; B; C;            % system matrices
        x;                  % current system state
        u;                  % current infusion rate (mg/mL*s)
        E;                  % current site effect (BIS)
        uMax;               % maximum allowed infusion rate (clinician-set)

        % pharmacokinetic parameters (drug)
        k10; k12; k13;      % transport coefficients
        k21; k31; ke0;     
        V1;                 % plasma deposition volume (L)

        % pharmacodynamic parameters (body)
        E0;                 % baseline site effect (BIS)
        Emax;               % desired upper bound of BIS
        ce50;               % site concentration at 50% BIS
        gamma;              % cooperativity constant in Hill formula

        % observer settings (extended kalman filter) 
        xEst; yEst;         % current estimate of concentrations and effect
        P; Q; R; S;         % estimatation error, process noise,
                            % sensor noise, and residual covariance
        H;                  % jacobian of measurement matrix (emax)
        K;                  % kalman gain, near optimal
        yError;             % residual of output effect
        errorBound;         % bound on residual error
        confidence;         % chi-squared confidence level
        chiSquared;         % chi-squared critical value

        % general controller settings
        ceMin;              % user-set site concentration lower bound
        Eop;                % operating (desired) site effect
        uMin;               % current minimum input from constraint
        w1; w2;             % weights for smoothing
        uOp;                % input required for the clinician-set
                            % operating point ceOp, assumed to be 50% BIS

        % exponential controller settings
        a1; a2;             % exponential cbf coefficients
        
        % graceful controller settings
        ceG;                % failsafe barrier
        delta;              % buffer defining new 
                            % safe set: ceG = ceMin + delta
        z;                  % zeta (damping ratio)
        w;                  % omega (natural frequency)

        % simulation history logs
        EHist;             % plant history
        xHist; uHist; 
        ceDotHist;
        tHist;              % discrete time vector for controller,
                            % observer, and discrete plant histories
        tHistContinuous;    % "continuous" time vector for continuous plant
        uMinHist;           % controller history
        xEstHist;           % observer history
        errorHist; 

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        % constructor 

        function sys = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x0, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w)

            % set initial values
            sys.P = P0; 

            % initialize plant parameters
            sys.k10 = k10; sys.k12 = k12; sys.k13 = k13; 
            sys.k21 = k21; sys.k31 = k31; sys.ke0 = ke0; 
            sys.V1 = V1; 
            sys.A = [-k10 - k12 - k13     k12     k13       0
                     k21                   -k21     0       0 
                     k31                    0      -k31     0 
                     ke0                    0       0      -ke0];
            sys.B = [1/V1; 0; 0; 0];
            sys.C = [0 0 0 1];
            sys.x = x0; sys.E0 = E0; sys.Emax = Emax; 
            sys.ce50 = ce50; sys.gamma = gamma;
            sys = sys.computeEffect(); % set initial E

            % initialize observer settings
            sys.Q = Q; sys.R = R; sys.confidence = confidence;
            sys.chiSquared = chi2inv(sys.confidence, 4);

            % initialize general controller settings   
            sys.uMin = 0; sys.uMax = uMax; 

            % initialize exponential controller settings
            sys.a1 = a1; sys.a2 = a2; 

            % initialize graceful controller settings
            sys.delta = delta; sys.z = z; sys.w = w; 

            % set bounds and op points
            sys.ceMin = ce50 * ((E0-Emax)/Emax)^(1/gamma);
            sys.ceG = sys.ceMin + sys.delta;
            ceOp = ce50 * ((E0-50)/50)^(1/gamma); % suppose clinician sets 50%
            sys.uOp = sys.V1 * sys.k10 * ceOp; 

            % initialize simulation history logs
            sys.tHist = zeros(1,1);
            sys.xHist = zeros(4, 1); sys.EHist = zeros(1,1); 
            sys.uHist = zeros(1,1); sys.uMinHist = zeros(1,1);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % plant function for ode45 integration
        
        function xDot = plant(sys, t, x, u)
        
            xDot = sys.A * x + sys.B * u; 

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extended kalman filter: a priori estimates

        function sys = predictEKF(sys, dt)

            % state prediction: euler discretization
            sys.xEst = sys.xEst + (sys.A * sys.xEst + sys.B * sys.u) * dt;

            % covariance prediction
            sys.P = sys.A * sys.P * sys.A' + sys.Q; 

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extended kalman filter: update step

        function sys = updateEKF(sys)

            % compute jacobian wrt a priori (1x4)
            sys.H = [sys.E0 * sys.gamma * ...
                (sys.ce50^(sys.gamma) * sys.xEst(4)^(sys.gamma-1)) / ...
                (sys.xEst(4)^(sys.gamma) + sys.ce50^(sys.gamma))^2 ...
                0 0 0];

            % compute residual e
            yTrue = sys.E; 
            sys.yEst = sys.E0 * (1 - ...
                (sys.xEst(4))^sys.gamma/ ...
                (sys.xEst(4)^sys.gamma + sys.ce50^sys.gamma));
            sys.yError= yTrue - sys.yEst;             

            % compute residual covariance S 
            sys.S = sys.H * sys.P * sys.H' + sys.R; 

            % compute (near-optimal) kalman gain
            sys.K = sys.P * sys.H' * (1/sys.S);

            % update a posteriori state estimate
            sys.xEst = sys.xEst + sys.K * sys.yError; 

            % update a posteriori covariance estimate
            sys.P = (eye(4) - sys.K * sys.H) * sys.P; 

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extended kalman filter: compute error bound based on confidence
        % ellipsoid

        function sys = computeErrorBoundEKF(sys)

            % find derivative of effect site concentration wrt effect
            deriv = -1 * sys.ce50 * (1/sys.gamma) * (sys.E0 / sys.E^2)...
                * ((sys.E0 - sys.E)/sys.E)^((1/sys.gamma) - 1); 

            % multiply derivative by current effect residual for
            % approximate effect site concentration residual
            sys.errorBound = deriv * sys.yError; 

            %fprintf("Approx. concentration residual: %.4f\n", sys.errorBound);
            

            % % compute eigenvalues of estimation error covariance
            % eigenvals = eigs(sys.P);
            % 
            % 
            % % compute semi-axis lengths
            % semiAxes = [sqrt(eigenvals(1) * sys.chiSquared); 
            %             sqrt(eigenvals(2) * sys.chiSquared); 
            %             sqrt(eigenvals(3) * sys.chiSquared); 
            %             sqrt(eigenvals(4) * sys.chiSquared)];
            % 
            % % find major axis 
            % sys.errorBound = max(semiAxes); 

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % exponential controller with linear constraint:
        % minimizes input itself + penalty for drastic changes,
        % subject to exponential constraint

        function sys = computeExponentialControl(sys, i, clamp)

            % compute current lower bound of input
            sys.uMin = (sys.V1/sys.ke0)*...
                (sys.ke0 * sys.k10 * sys.x(1) - ...
                sys.ke0 * sys.k12 * (sys.x(2) - sys.x(1)) - ...
                sys.ke0 * sys.k13 * (sys.x(3) - sys.x(1)) - ...
                (sys.a1 + sys.a2) * (sys.ke0 * (sys.x(1) - sys.x(4))) - ...
                sys.a1 * sys.a2 * (sys.x(4) - sys.ceMin));

            if (sys.uMin > sys.uMax)
                disp('   Over actuation limit!');
                if clamp == true
                    sys.u = sys.uMax; 
                end
            elseif (sys.uMin > sys.uOp)
                sys.u = sys.uMin;
            else
                sys.u = sys.uOp;
            end

        end

        function uMin = computeExponentialMin(sys)
            uMin = (sys.V1/sys.ke0)*...
                (sys.ke0 * sys.k10 * sys.x(1) - ...
                sys.ke0 * sys.k12 * (sys.x(2) - sys.x(1)) - ...
                sys.ke0 * sys.k13 * (sys.x(3) - sys.x(1)) - ...
                (sys.a1 + sys.a2) * (sys.ke0 * (sys.x(1) - sys.x(4))) - ...
                sys.a1 * sys.a2 * (sys.x(4) - sys.ceMin));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % graceful controller with nonlinear constraint:
        % minimizes squared difference of input from uOp, 
        % subject to nonlinear graceful constraint. 
        % utilizes error bound on EKF for robustness. 

        function sys = computeGracefulControl(sys, clamp)

            hg = (sys.x(4)-sys.ceG)/(sys.ceG - sys.ceMin);
            hgDot = sys.ke0 * (sys.x(1)-sys.x(4))/(sys.ceG-sys.ceMin);

            % compute current lower bound of input, based on nonlinear
            % mass-spring damper constraint 

            sys.uMin = ...
                sys.x(1) * sys.V1 * (sys.k10 + sys.k12 + sys.k13 ...
                    + sys.ke0)...
              - sys.x(2) * sys.V1 * sys.k12...
              - sys.x(3) * sys.V1 * sys.k13...
              - sys.x(4) * sys.V1 * sys.ke0...
              - sys.V1 * ((sys.ceG-sys.ceMin)/sys.ke0) * ...
                ((2 * sys.w * sys.z * hgDot) + sys.w^2 * (hg / (hg + 1)) ...
                + abs(sys.errorBound));

            % compute control input with switching control policy

            if (sys.uMin > sys.uMax)
                fprintf('   Over actuation limit!');
                if clamp == true
                    sys.u = sys.uMax;
                    fprintf(' Clamped.')
                end
                fprintf('\n');
            elseif (sys.uMin > sys.uOp)
                sys.u = sys.uMin; 
            else
                sys.u = sys.uOp;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % update the e-max site effect as a function of current effect site
        % concentration

        function sys = computeEffect(sys)

            % compute current site effect
            sys.E = sys.E0 * (1 - ...
                (sys.x(4))^sys.gamma/ ...
                (sys.x(4)^sys.gamma + sys.ce50^sys.gamma));

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % advance continuous simulation by one step

        function sys = step(sys, i, t0, dt, graceful, clamp)

            % observer: prediction and update
            sys = sys.predictEKF(dt);
            sys = sys.updateEKF;
            sys = sys.computeErrorBoundEKF;

            % controller: compute and apply control input
            if (graceful == true)
                sys = sys.computeGracefulControl(clamp); 
            else
                sys = sys.computeExponentialControl(i, clamp);
            end

            % plant: advance (continuous)
            integrand = @(t, x) sys.plant(t, x, sys.u);
            [t, xt] = ode45(integrand, [t0, t0 + dt], sys.x);
            sys.x = xt(end, :);
            sys.tHistContinuous = [sys.tHistContinuous, t'];
            sys = sys.computeEffect();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % simulate continuous plant and discrete controller

        function [sys, tHist, tHistContinuous, xHist, xEstHist,...
                errorHist, EHist, uHist, ...
                uMinHist, ceDotHist] = simulate(sys, ...
                graceful, clamp, x0, xEst0, P0, tStart, tEnd, nTimes)

            % compute discrete time interval
            dt = (tEnd - tStart) / nTimes;

            % reset history vectors
            sys.tHist = zeros(1, nTimes); sys.xEstHist = zeros(4, nTimes);
            sys.errorHist = zeros(1, nTimes);
            sys.xHist = zeros(4, nTimes); sys.EHist = zeros(1, nTimes); 
            sys.uHist = zeros(1, nTimes); sys.uMinHist = zeros(1, nTimes);
            sys.ceDotHist = zeros(1, nTimes);

            % set and log initial conditions
            sys.x = x0; sys.xEst = xEst0; sys.P = P0; 
            sys.u = 0;
            sys.tHist(1) = tStart;
            sys.xHist(:,1) = x0; 
            sys.ceDotHist = sys.ke0 * (sys.x(1) - sys.x(4));

            t0 = tStart;
            
            for i=2:nTimes

                % time step
                sys = sys.step(i, t0, dt, graceful, clamp);

                % log
                sys.tHist(i) = sys.tHist(i-1) + dt;
                sys.xHist(:,i) = sys.x; sys.xEstHist(:,i) = sys.xEst; 
                sys.errorHist(i) = sys.yError; sys.EHist(i) = sys.E;
                sys.uHist(i) = sys.u; sys.uMinHist(i) = sys.uMin;
                sys.ceDotHist(i) = sys.ke0 * (sys.x(1) - sys.x(4));

                % update t0
                t0 = t0 + dt; 

            end

            tHist = sys.tHist; tHistContinuous = sys.tHistContinuous;
            xHist = sys.xHist; xEstHist = sys.xEstHist; EHist = sys.EHist; 
            uHist = sys.uHist; uMinHist = sys.uMinHist;
            ceDotHist = sys.ceDotHist; errorHist = sys.errorHist; 

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % plot the characteristic e-max vs. effect site concentration
        % curve of this system (i.e. this individual patient)

        function plotCharacteristicCurve(sys)

            sys.setPlotSettings();     

            % plot characteristic effect curve (not simulated)
            ce = linspace(0, 100);
            curve = zeros(100);
            for i=1:100
                curve(i) = sys.E0 * (1 - ...
                ce(i)^sys.gamma /(ce(i)^sys.gamma + sys.ce50^sys.gamma));
            end
            figure('Color', [1 1 1]);
            plot(ce, ...
                curve);
            xlim([0 6]);
            title('Patient characteristic propofol effect curve');
            xlabel('Effect site propofol concentration (ug/mL)')
            ylabel('Propofol effect (BIS)')
            ax = gca;
            ax.TitleFontSizeMultiplier = 1.5;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % makes 2 plots: (1) estimated vs. true states over time; 
        % (2) state and output estimation residuals over time

        function plotEstimation(sys)

            sys.setPlotSettings(); 
    
            % plot estimated vs. true states
            figure('Color', [1 1 1]);
            hold on; 
            plot(sys.tHist, sys.xEstHist(1,:), '--'); 
            plot(sys.tHist, sys.xEstHist(2,:), '--'); 
            plot(sys.tHist, sys.xEstHist(3,:), '--'); 
            plot(sys.tHist, sys.xEstHist(4,:), '--'); 
            plot(sys.tHist, sys.xHist(1,:), '-'); 
            plot(sys.tHist, sys.xHist(2,:), '-'); 
            plot(sys.tHist, sys.xHist(3,:), '-'); 
            plot(sys.tHist, sys.xHist(4,:), '-'); 
            legend('Plasma concentration', ...
                    'Fast peripheral concentration', ...
                    'Slow peripheral concentration', ...
                    'Effect site concentration');
            title('EKF Estimated vs. True States')
            xlabel('Time [min]');
            ylabel('Concentration [ug mL⁻¹]')


            % plot error over time
            figure('Color', [1 1 1]);
            hold on; 
            plot(sys.tHist, sys.errorHist); 
            plot(sys.tHist, sys.xEstHist(1,:) - sys.xHist(1,:));
            plot(sys.tHist, sys.xEstHist(2,:) - sys.xHist(2,:));
            plot(sys.tHist, sys.xEstHist(3,:) - sys.xHist(3,:));
            plot(sys.tHist, sys.xEstHist(4,:) - sys.xHist(4,:));
            legend('BIS error', 'Plasma concentration error', ...
                'Fast peripheral concentration error', ...
                'Slow peripheral concentration error', ...
                'Effect site concentration error');
            title('EKF Residuals');
            xlabel('Time [min]');
            ylabel('Error in concentration [ug mL⁻¹]')
        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % plot phase plane showing xHist trajectory and initial actuation
        % boundary

        function plotActuationLimit(sys)

            sys.setPlotSettings;

            figure();
            plot(ceHist, ceDostHist); 
            

        end

    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static)

        % cosmetic settings

        function setPlotSettings()

            set(groot, 'defaultTextInterpreter', 'tex');
            set(groot, 'defaultLineLineWidth', 2);
            set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2);  
            set(groot, 'defaultAxesFontSize', 12);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end