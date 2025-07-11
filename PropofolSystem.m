classdef PropofolSystem

    % an automated intravenous propofol infusion control system
    % using a three-compartment mamillary pkpd model. 
    % implements an extended kalman filter and barrier function
    % controllers. 

    properties

        % plant
        A; B; C;            % system matrices
        x; x0;              % current and initial system states
        u;                  % current infusion rate [ug mL⁻¹ min⁻¹]
        BIS;                % current BIS [%]

        % pharmacokinetic parameters
        k10; k12; k13;      % transport coefficients
        k21; k31; ke0;     
        V1;                 % plasma deposition volume (L)

        % pharmacodynamic parameters
        BISdes;             % desired BIS [%]
        BISmax;             % desired maximum of BIS [%]
        BIS0;               % baseline drugless BIS [%]
        ce50;               % site concentration at 50% BIS [ug mL⁻¹]
        gamma;              % Hill cooperativity constant []

        % observer settings (extended kalman filter) 
        xEst; BISest;       % current estimate of concentrations and BIS
        P; Q; R; S;         % estimatation error, process noise,
                            % sensor noise, and residual covariance
        H;                  % jacobian of measurement matrix (BISmax)
        K;                  % kalman gain, near optimal
        errorBound;         % bound on residual error
        confidence;         % chi-squared confidence level
        chiSquared;         % chi-squared critical value

        % general controller settings
        ceMin;              % set concentration lower bound
        ceDes;              % desired concentration
        uMin;               % current minimum input from constraint
        uOp;                % input required for clinician-set BISdes
        uMax;               % maximum allowed infusion rate (clinician-set)

        % exponential controller settings
        a1; a2;             % exponential cbf coefficients
        
        % graceful controller settings
        ceG;                % failsafe barrier
        z;                  % zeta (damping ratio)
        w;                  % omega (natural frequency)

        % simulation history logs
        BISHist;            % plant history
        xHist; uHist; 
        ceDotHist;
        tHist;              % discrete time vector for controller,
                            % observer, and discrete plant histories
        tHistContinuous;    % "continuous" time vector for continuous plant
        uMinHist;           % controller history
        xEstHist;           % observer history
        BISEstHist; 

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        % constructor 

        function sys = PropofolSystem(...
                x0, k10, k12, k13, k21, k31, ke0, V1, ...  
                BIS0, BISdes, BISmax, ce50, gamma, ...  
                xEst0, P0, Q, R, ...
                uMax, a1, a2, delta, z, w)

            % initialize pharmacokinetic plant parameters
            sys.x = x0; sys.x0 = x0; 
            sys.k10 = k10; sys.k12 = k12; sys.k13 = k13; 
            sys.k21 = k21; sys.k31 = k31; sys.ke0 = ke0; 
            sys.V1 = V1; 
            sys.A = [-k10 - k12 - k13     k12     k13       0
                     k21                 -k21     0         0 
                     k31                    0    -k31       0 
                     ke0                    0      0       -ke0];
            sys.B = [1/V1; 0; 0; 0];
            sys.C = [0 0 0 1];

            % initialize pharmacodynamic parameters
            sys.BIS0 = BIS0; sys.ce50 = ce50; sys.gamma = gamma;
            sys = sys.computeEffect(); 
            
            % initialize general controller settings   
            sys.uMax = uMax; sys.BISdes = BISdes; sys.BISmax = BISmax;

            % use clinician-set BIS operating point for ce operating point
            sys.ceDes = ce50 * ((BIS0 - BISdes)/BISdes)^(1/gamma);
            sys.uOp = sys.V1 * sys.k10 * sys.ceDes;
            
            % use clinician-set BIS maximum for ce minimum
            sys.ceMin = ce50 * ((BIS0 - BISmax)/BISmax)^(1/gamma);

            % initialize exponential controller settings
            sys.a1 = a1; sys.a2 = a2; 

            % initialize graceful controller settings
            sys.z = z; sys.w = w;
            sys.ceG = sys.ceMin + delta;

            % initialize observer settings
            sys.xEst = xEst0; sys.P = P0; 
            sys.Q = Q; sys.R = R;

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
            sys.H = [0 0 0 ...
                sys.BIS0 * sys.gamma * ...
                (sys.ce50^(sys.gamma) * sys.xEst(4)^(sys.gamma-1)) / ...
                (sys.xEst(4)^(sys.gamma) + sys.ce50^(sys.gamma))^2];

            % compute output residual
            sys.BISest = sys.BIS0 * (1 - ...
                (sys.xEst(4))^sys.gamma/ ...
                (sys.xEst(4)^sys.gamma + sys.ce50^sys.gamma));
            BISerror = sys.BIS - sys.BISest;  

            % compute residual covariance S 
            sys.S = sys.H * sys.P * sys.H' + sys.R; 

            % compute (near-optimal) kalman gain
            sys.K = sys.P * sys.H' * (1/sys.S);

            % update a posteriori state estimate
            sys.xEst = sys.xEst + sys.K * BISerror; 

            % update a posteriori covariance estimate
            sys.P = (eye(4) - sys.K * sys.H) * sys.P; 

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extended kalman filter: compute error bound based on confidence
        % ellipsoid

        function sys = computeErrorBoundEKF(sys)

            % find derivative of effect site concentration wrt effect
            deriv = sys.BIS0 * sys.gamma * ...
                (sys.ce50^(sys.gamma) * sys.xEst(4)^(sys.gamma-1)) / ...
                (sys.xEst(4)^(sys.gamma) + sys.ce50^(sys.gamma))^2 ;

            % multiply derivative by current effect residual for
            % approximate effect site concentration residual
            sys.errorBound = sys.BIS / deriv; 

            %fprintf("Approx. concentration residual: %.4f\n", sys.BISrrorBound);
            

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
            % sys.BISrrorBound = max(semiAxes); 

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % exponential controller with linear constraint:
        % minimizes input itself + penalty for drastic changes,
        % subject to exponential constraint

        function sys = computeExponentialControl(sys, clamp, observer)

            % compute current lower bound of input
            if observer == true
                sys.uMin = (sys.V1/sys.ke0)*...
                (sys.ke0 * sys.k10 * sys.xEst(1) - ...
                sys.ke0 * sys.k12 * (sys.xEst(2) - sys.xEst(1)) - ...
                sys.ke0 * sys.k13 * (sys.xEst(3) - sys.xEst(1)) - ...
                (sys.a1 + sys.a2) * (sys.ke0 * (sys.xEst(1) - ...
                sys.xEst(4))) - sys.a1 * sys.a2 * (sys.xEst(4) - sys.ceMin));
            else
                sys.uMin = (sys.V1/sys.ke0)*...
                (sys.ke0 * sys.k10 * sys.x(1) - ...
                sys.ke0 * sys.k12 * (sys.x(2) - sys.x(1)) - ...
                sys.ke0 * sys.k13 * (sys.x(3) - sys.x(1)) - ...
                (sys.a1 + sys.a2) * (sys.ke0 * (sys.x(1) - sys.x(4))) - ...
                sys.a1 * sys.a2 * (sys.x(4) - sys.ceMin));
            end
            
            % kkt solution
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % graceful controller with nonlinear constraint:
        % minimizes squared difference of input from uOp, 
        % subject to nonlinear graceful constraint. 
        % utilizes error bound on EKF for robustness. 

        function sys = computeGracefulControl(sys, clamp, observer)

            % compute current lower bound of input, based on nonlinear
            % mass-spring damper constraint

            hg = (sys.x(4)-sys.ceG)/(sys.ceG-sys.ceMin);
            hgDot = sys.ke0 * (sys.x(1) - sys.x(4)) / (sys.ceG-sys.ceMin);

            sys.uMin = ...
                sys.x(1) * sys.V1 * (sys.k10 + sys.k12 + sys.k13 ...
                    + sys.ke0)...
              - sys.x(2) * sys.V1 * sys.k12...
              - sys.x(3) * sys.V1 * sys.k13...
              - sys.x(4) * sys.V1 * sys.ke0...
              - sys.V1 * ((sys.ceG-sys.ceMin)/sys.ke0) * ...
                ((2 * sys.w * sys.z * hgDot) + sys.w^2 * (hg / (hg + 1)));

            % compute control input with switching control policy

            if (sys.uMin > sys.uMax)
                fprintf('   Over actuation limit!');
                if clamp == true
                    sys.u = sys.uMax;
                    fprintf(' Clamped.')
                else 
                    sys.u = sys.uMin; 
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
            sys.BIS = sys.BIS0 * (1 - ...
                (sys.x(4))^sys.gamma/ ...
                (sys.x(4)^sys.gamma + sys.ce50^sys.gamma));

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % advance continuous simulation by one step

        function sys = step(sys, t0, dt, graceful, clamp, observer)

            % observer: prediction and update
            sys = sys.predictEKF(dt);
            sys = sys.updateEKF;
            sys = sys.computeErrorBoundEKF;

            % controller: compute and apply control input
            if (graceful == true)
                sys = sys.computeGracefulControl(clamp, observer); 
            else
                sys = sys.computeExponentialControl(clamp, observer);
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

        function [sys, tHist, tHistContinuous, ...
                xHist, xEstHist,...
                uHist, uMinHist, ceDotHist,...
                BISHist, BISEstHist] = simulate(sys, ...
                graceful, clamp, observer, ...
                tStart, tEnd, nTimes)

            % compute discrete time interval
            dt = (tEnd - tStart) / nTimes;

            % reset history vectors
            sys.tHist = zeros(1, nTimes); sys.xEstHist = zeros(4, nTimes);
            sys.BISEstHist = zeros(1, nTimes);
            sys.xHist = zeros(4, nTimes); sys.xEstHist = zeros(4, nTimes);
            sys.BISHist = zeros(1, nTimes); sys.BISEstHist = zeros(1, nTimes);
            sys.uHist = zeros(1, nTimes); sys.uMinHist = zeros(1, nTimes);
            sys.ceDotHist = zeros(1, nTimes);

            % first logs
            sys.x = sys.x0; sys.u = 0;
            sys.tHist(1) = tStart;
            sys.xHist(:,1) = sys.x; 
            sys.ceDotHist(1) = dot(sys.A(4,:), sys.x);

            t0 = tStart;
            
            for i=2:nTimes

                % time step
                sys = sys.step(t0, dt, graceful, clamp, observer);

                % log
                sys.tHist(i) = sys.tHist(i-1) + dt;
                sys.xHist(:,i) = sys.x; sys.xEstHist(:,i) = sys.xEst; 
                sys.BISHist(i) = sys.BIS; sys.BISEstHist(i) = sys.BISest; 
                sys.uHist(i) = sys.u; sys.uMinHist(i) = sys.uMin;
                sys.ceDotHist(i) = sys.ke0 * (sys.x(1) - sys.x(4));

                % update t0
                t0 = t0 + dt; 

            end

            tHist = sys.tHist; tHistContinuous = sys.tHistContinuous;
            xHist = sys.xHist; xEstHist = sys.xEstHist;
            uHist = sys.uHist; uMinHist = sys.uMinHist;
            BISEstHist = sys.BISEstHist; BISHist = sys.BISHist; 
            ceDotHist = sys.ceDotHist; 

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
                curve(i) = sys.BIS0 * (1 - ...
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
            legend('Plasma concentration estimate', ...
                    'Fast peripheral concentration estimate', ...
                    'Slow peripheral concentration estimate', ...
                    'Effect site concentration estimate', ...
                    'Plasma concentration true state',...
                    'Fast peripheral concentration true state',...
                    'Slow peripheral concentration true state');
            title('EKF Estimated vs. True States')
            xlabel('Time [min]');
            ylabel('Concentration [ug mL⁻¹]')

            % plot error over time
            figure('Color', [1 1 1]);
            hold on; 
            plot(sys.tHist, sys.BISEstHist - sys.BISHist); 
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
            set(groot, 'defaultAxesFontSize', 15);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end