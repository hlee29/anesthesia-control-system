classdef PropofolSystem

    % A patient modeled as a PK/PD system for surgical propofol infusion.
    % Methods: 
    %       PropofolSystem() constructor
    %       xDot() for ode15s integrand
    %       computeGracefulControl()
    %       getConservativeGracefulBounds()
    %       estimateInterval()
    %       step()
    %       simulate()
    %       plotEstimation()
    %       plotCharacteristicCurve()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        
        % Plot color
        color; 

        % Plant
        A; B; C;                            % System matrices
        x; x0;                              % Current, initial states
        u;                                  % Current infusion rate 
        BIS;                                % Current BIS
        k10; k12; k13; k21;                 % PK parameters
        k31; ke0; V1; 
        BIS0; ce50; gamma;                  % PD parameters

        % Controller constraints
        uCap;                               % Set actuation limit
        uOp;                                % Desired input
        uMin; uMax;                         % Computed input bounds
        BISdes; BISmin; BISmax;             % Desired BIS bounds
        ceDes; ceMin; ceMax;                % Desired effect site conc.

        % Controller parameters
        ceG1; ceG2;                         % Lower and upper barriers
        zeta;                               % CBF damping ratio
        omega;                              % CBF natural frequency
        wGrc1; wGrc2;                       % Min, max slack pareto weights
        ceMinGrc;                           % Lowermost barrier

        % Interval observer parameters
        xLoEst0; xHiEst0;                   % Initial estimates
        xLoEst; xHiEst;                     % Current estimates
        zLoEst0; zHiEst0;                   % Current estimates (diagonal basis)
        zLoEst; zHiEst;
        dLo, dHi;                           % Current disturbances
        V;                                  % Current noise bounds
        L;                                  % Gain
        T; D; Tinv;                         % Diagonalization 

        % Simulation history data
        xHist; uHist; BISHist;              % True states, input, output
        ceDotHist;                          % Derivative of x4
        tHist; tHistContinuous;             % Time vectors
        xLoEstHist; xHiEstHist;             % Estimated state intervals
        vHist; dHist;                       % Disturbance and noise

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        % Constructor
        function sys = PropofolSystem(color, x0, ...
                k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, ce50, gamma, ...
                BISdes, BISmin, BISmax, ...
                uCap, a1, a2, wExp, ...
                ceMinGrc, zeta, omega, wGrc1, wGrc2, ...
                xLoEst0, xHiEst0, d, V)
            sys.color = color; sys.x0 = x0; 
            sys.k10 = k10; sys.k12 = k12; 
            sys.k13 = k13; sys.k21 = k21;  
            sys.k31 = k31; sys.ke0 = ke0; sys.V1 = V1; 
            sys.BIS0 = BIS0; sys.ce50 = ce50; sys.gamma = gamma; 
            sys.BISdes = BISdes; sys.BISmin = BISmin; sys.BISmax = BISmax; 
            sys.ceDes = sys.invertBIS(BISdes);
            sys.ceMin = sys.invertBIS(BISmax); 
            sys.ceMax = sys.invertBIS(BISmin); 
            sys.uCap = uCap; sys.a1 = a1; sys.a2 = a2; sys.wExp = wExp; 
            sys.ceG1 = sys.ceDes; 
            sys.ceG2 = sys.invertBIS(BISdes-5); sys.wGrc1 = wGrc1; 
            sys.wGrc2 = wGrc2; sys.ceMinGrc = ceMinGrc;
            sys.zeta = zeta; sys.omega = omega;
            sys.xLoEst0 = xLoEst0; sys.xHiEst0 = xHiEst0;
            
            sys.dLo = -d; sys.dHi = d; 
            sys.V = V; 

            % Construct system matrices
            sys.A = [-k10-k12-k13   k12     k13     0
                      k21          -k21     0       0
                      k31           0      -k31     0
                      ke0           0       0       -ke0];
            sys.B = [1/V1; 0; 0; 0;]; 
            sys.C = [0 0 0 1]; 
            sys.L = [0.9724; 2.3618; 0.1399; 4.8826];
            sys.T = [0.1742    0.9993   -0.3699   -0.2252
                    0.4232    0.0176   -0.9290   -0.0000
                    0.0248    0.0010   -0.0003   -0.9743
                    0.8888    0.0342   -0.0109   -0.0065]; 

            % Observer basis
            sys.Tinv = inv(sys.T);
            sys.D = sys.Tinv * (sys.A - sys.L * sys.C) * sys.T; 
            sys.zLoEst0 = sys.Tinv * sys.xLoEst0;
            sys.zHiEst0 = sys.Tinv * sys.xHiEst0;
            sys.zLoEst = sys.zLoEst0;
            sys.zHiEst = sys.zHiEst0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Plant function for ode45 integration
        function xDot = xDot(sys, t, x, u)
            xDot = sys.A*x + sys.B*u; 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = computeGracefulControl(sys, useObserver)
            if ~useObserver
                hgMin = (sys.x(4)-sys.ceG1) / (sys.ceG1-sys.ceMinGrc);
                hgDotMin = sys.ke0 * (sys.x(1)-sys.x(4)) / ...
                    (sys.ceG1-sys.ceMinGrc);
                sys.uMin = sys.x(1) * sys.V1 * (sys.k10 + sys.k12 + ...
                    sys.k13 + sys.ke0)...
                      - sys.x(2) * sys.V1 * sys.k12...
                      - sys.x(3) * sys.V1 * sys.k13...
                      - sys.x(4) * sys.V1 * sys.ke0...
                      - sys.V1 * ((sys.ceG1-sys.ceMinGrc)/sys.ke0) * ...
                        ((2 * sys.omega * sys.zeta * hgDotMin) + ...
                        sys.omega^2 * (hgMin / (hgMin+1)));

                hgMax = (sys.ceG2-sys.x(4))/(sys.ceMax-sys.ceG2); 
                hgDotMax = sys.ke0 * (sys.x(4)-sys.x(1)) / ...
                    (sys.ceMax-sys.ceG2); 
                sys.uMax = sys.V1 * sys.ke0 * (sys.x(1)-sys.x(4)) - ...
                    sys.V1 * sys.k12 * (sys.x(2)-sys.x(1)) - ...
                    sys.V1 * sys.k13 * (sys.x(3)-sys.x(1)) + ...
                    sys.V1 * sys.k10 * sys.x(1) + ...
                    sys.V1 * ((sys.ceMax-sys.ceG2)/sys.ke0) * ...
                    (2*sys.zeta*sys.omega*hgDotMax + ...
                    sys.omega^2 * (hgMax / (hgMax+1)));
                sys.uMax = max(0, sys.uMax);
            else 
                [sys.uMin, sys.uMax] = sys.getConservativeGracefulBounds();
            end
                        
            % Solve quadprog [u, s1, s2]
            H = [1 0 0; 0 sys.wGrc1 0; 0 0 sys.wGrc2]; f = [0; 0; 0];  
            A = [1 1 0; -1 0 -1]; b = [sys.uMax; -1 * sys.uMin];
            lb = [0, 0, 0]; ub = [sys.uCap, inf, inf];
            Aeq = []; beq = []; x0 = []; 
            options = optimoptions('quadprog', 'Display', 'off');
            solution = quadprog(H, f, A, b, Aeq, beq, lb, ub, x0, options); 
            sys.u = solution(1); 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute conservative safety bounds on the graceful control input
        % given state interval estimation

        function [uMinConservative, uMaxConservative] = ...
        getConservativeGracefulBounds(sys)

            % Maximize uMin given xEst in [xLoEst, xHiEst]
            uMinConservative = ...
                (sys.k10 + sys.k12 + sys.ke0 + ...
                    (2*sys.V1*sys.omega*sys.zeta)) * sys.xHiEst(1) ...
                - sys.k12 * sys.xLoEst(2) ...
                - sys.k13 * sys.xLoEst(3) ...
                - (sys.ke0 + ...
                    (2*sys.V1*sys.omega*sys.zeta)) * sys.xLoEst(4) ...
                - sys.V1 * ((sys.ceG1-sys.ceMinGrc)/sys.ke0) * ...
                    sys.omega^2 * ((sys.xLoEst(4)-sys.ceG1)/ ...
                        (sys.xLoEst(4)-sys.ceMinGrc)); 

            % Minimize uMax given xEst in [xLoEst, xHiEst]
            uMaxConservative = ...
                (sys.ke0 + sys.k12 + sys.k13 + sys.k10 - ...
                    (2*sys.V1*sys.omega*sys.zeta)) * sys.xHiEst(1) ...
               - sys.k12 * sys.xHiEst(2) ...
               - sys.k13 * sys.xHiEst(3) ... 
               + 2*sys.V1*sys.omega*sys.zeta * sys.xLoEst(4) ...
               + sys.V1 * ((sys.ceMax-sys.ceG2)/sys.ke0) * ...
                    sys.omega^2 * ((sys.ceG2-sys.xHiEst(4))/ ...
                        (sys.ceMax-sys.xHiEst(4)));
            uMaxConservative = max(0, uMaxConservative);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute the bispectral index as a Hill function of the 
        % effect site concentration 

        function BIS = computeBIS(sys, ce)
            BIS = sys.BIS0 * (1 - ...
                (ce^sys.gamma)/ ...
                (ce^sys.gamma + sys.ce50^sys.gamma));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute the effect site concentration as an inverted Hill
        % function of the bispectral index

        function ceMeasured = invertBIS(sys, BIS)
            ceMeasured = sys.ce50 * ...
                ((sys.BIS0 - BIS) / BIS)^(1/sys.gamma);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Update interval estimator's estimates

        function sys = estimateInterval(sys, dt, i)

            % Compute z-bounds
            ceMeasured = sys.invertBIS(sys.BIS); 
            Bz = sys.Tinv * sys.B; 
            Lz = sys.Tinv * sys.L; 
            distBoundZ = abs(sys.Tinv) * (ones(4,1) * sys.dHi); 
            noiseBoundZ = abs(Lz) * sys.V; 
            sys.zLoEst = sys.zLoEst + ...
                dt * (real(sys.D) * sys.zLoEst + Bz * sys.u + ...
                    Lz * ceMeasured - noiseBoundZ - distBoundZ);
            sys.zHiEst = sys.zHiEst + ...
                dt * (real(sys.D) * sys.zHiEst + Bz * sys.u + ...
                    Lz * ceMeasured + noiseBoundZ + distBoundZ); 

            % Recover x-bounds, minding positivity of concentration
            Tpos = max(sys.T, 0); Tneg = min(sys.T, 0);
            sys.xHiEst = Tpos * sys.zHiEst + Tneg * sys.zLoEst;
            sys.xLoEst = Tpos * sys.zLoEst + Tneg * sys.zHiEst;
            sys.xLoEst = max(0, sys.xLoEst);  
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Advance a step in the simulation

        function sys = step(sys, t0, dt, graceful, runObserver, ...
                useObserver, i)

            % Advance observer
            if runObserver
                sys = sys.estimateInterval(dt, i); 
            end

            % Advance controller
            if graceful
                sys = sys.computeGracefulControl(useObserver); 
            else 
                sys = sys.computeExponentialControl(); 
            end

            % Advance plant
            integrand = @(t,x) sys.xDot(t, x, sys.u); 
            [t, xt] = ode15s(integrand, [t0, t0+dt], sys.x); 
            sys.x = xt(end,:); 
            sys.tHistContinuous = [sys.tHistContinuous, t']; 
            sys.BIS = sys.computeBIS(sys.x(4)); 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Simulate the system

        function sys = simulate(sys, graceful, runObserver, ...
                useObserver, tStart, tEnd, nTimes)

            % Compute discrete time interval
            dt = (tEnd-tStart) / nTimes;

            % Reset history vectors
            sys.tHist = zeros(1, nTimes); 
            sys.xHist = zeros(4, nTimes); sys.ceDotHist = zeros(1, nTimes); 
            sys.xLoEst = sys.xLoEst0; sys.xHiEst = sys.xHiEst0;
            sys.xLoEstHist = zeros(4, nTimes); 
            sys.xLoEstHist(:,1) = sys.xLoEst;
            sys.xHiEstHist = zeros(4, nTimes); 
            sys.xHiEstHist(:,1) = sys.xHiEst;
            sys.uHist = zeros(1, nTimes); 
            sys.BISHist = zeros(1, nTimes); 
            sys.dHist = zeros(8, nTimes);

            % First logs
            sys.x = sys.x0; sys.xHist(:,1) = sys.x; 
            sys.ceDotHist(1) = dot(sys.A(4,:), sys.x); 
            sys.xLoEstHist(:,1) = sys.xLoEst0; 
            sys.xHiEstHist(:,1) = sys.xHiEst0; 
            t0 = tStart; sys.u = 0; sys.BIS = sys.computeBIS(sys.x(4)); 
            sys.BISHist(1) = sys.BIS; 
            for i = 2:nTimes
                sys = sys.step(t0, dt, graceful, runObserver, ...
                    useObserver, i); 
                sys.tHist(i) = sys.tHist(i-1) + dt; 
                sys.xHist(:,i) = sys.x;
                sys.ceDotHist(i) = dot(sys.A(4,:), sys.x);
                sys.xLoEstHist(:,i) = sys.xLoEst; 
                sys.xHiEstHist(:,i) = sys.xHiEst; 
                sys.BISHist(i) = sys.BIS; sys.uHist(i) = sys.u; 
                sys.dHist(1:4,i) = sys.dLo; sys.dHist(5:8,i) = sys.dHi; 
                t0 = t0 + dt; 
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Plot the patient's characteristic curve of bispectral index
        % with respect to effect site concentration

        function plotCharacteristicCurve(sys)
            sys.setPlotSettings(); 
            ce = linspace(0, 100); 
            curve = zeros(100); 
            for i = 1:100
                curve(i) = sys.computeBIS(ce(i)); 
            end
            figure('Color', [1 1 1]); 
            plot(ce, curve); 
            xlim([0, 6]);  
            title('Patient characteristic propofol effect curve'); 
            xlabel('Effect site propofol concentration [μg mL⁻¹]'); 
            ylabel('Bispectral index [%]'); 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Plot the patient's true vs. estimated trajectories over time

        function plotEstimation(sys, num)
            rgb = sscanf(sys.color(2:end),'%2x%2x%2x',[1 3])/255;
            sys.setPlotSettings(); 
            figure('Color', [1 1 1]); 
            hold on
            subplot(2, 2, 1)
                hold on
                fill([sys.tHist, fliplr(sys.tHist)], ...
                    [sys.xHiEstHist(1,:), fliplr(sys.xLoEstHist(1,:))], ...
                    rgb, 'FaceAlpha', 0.4, 'EdgeColor','none');
                plot(sys.tHist, sys.xHist(1,:), '-', 'Color', rgb); 
                title('Plasma concentration', 'FontSize', 12);
                xlabel('Time [min]', 'FontSize', 12); 
                ylabel('Concentration [μg mL⁻¹]', 'FontSize', 12);
            subplot(2, 2, 2)
                hold on
                fill([sys.tHist, fliplr(sys.tHist)], ...
                     [sys.xHiEstHist(2,:), fliplr(sys.xLoEstHist(2,:))], ...
                     rgb, 'FaceAlpha', 0.4, 'EdgeColor','none');
                plot(sys.tHist, sys.xHist(2,:), '-', 'Color', rgb); 
                title('Fast peripheral concentration', 'FontSize', 12);
                xlabel('Time [min]', 'FontSize', 12); 
                ylabel('Concentration [μg mL⁻¹]', 'FontSize', 12);
            subplot(2, 2, 3)
                hold on
                fill([sys.tHist, fliplr(sys.tHist)], ...
                    [sys.xHiEstHist(3,:), fliplr(sys.xLoEstHist(3,:))], ...
                    rgb, 'FaceAlpha', 0.4, 'EdgeColor','none');
                plot(sys.tHist, sys.xHist(3,:), '-', 'Color', rgb);
                title('Slow peripheral concentration', 'FontSize', 12);
                xlabel('Time [min]', 'FontSize', 12); 
                ylabel('Concentration [μg mL⁻¹]', 'FontSize', 12);
            subplot(2, 2, 4)
                hold on
                fill([sys.tHist, fliplr(sys.tHist)], ...
                    [sys.xHiEstHist(4,:), fliplr(sys.xLoEstHist(4,:))], ...
                     rgb, 'FaceAlpha', 0.4, 'EdgeColor','none');
                plot(sys.tHist, sys.xHist(4,:), '-', 'Color', rgb);
                yline(sys.ceMin, 'k--', 'safety min', 'FontSize', 10);
                yline(sys.ceMax, 'k--', 'safety max', 'FontSize', 10);
                yline(sys.ceDes, 'k--', 'desired', 'FontSize', 10);
                title('Effect site concentration', 'FontSize', 12);
                xlabel('Time [min]', 'FontSize', 12);
                ylabel('Concentration [μg mL⁻¹]', 'FontSize', 12);
            sgtitle(['True states vs. Estimated intervals for Patient ' num]);            
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static)
        function setPlotSettings()
            set(groot, 'defaultTextInterpreter', 'tex'); 
            set(groot, 'defaultLineLineWidth', 2); 
            set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2); 
            set(groot, 'defaultAxesFontSize', 15); 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
