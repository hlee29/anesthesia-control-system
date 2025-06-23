classdef PropofolSystem

    % an automated IV propofol infusion control system

    properties

        % linear plant: infusion mass/vol input, site concentration output
        A; B; C;            % system matrices
        x;                  % current system state
        u;                  % current input (mg/mL*s)
        E;                  % current site effect (BIS)

        % optimal controller
        ceMin;              % site concentration lower bound
        Eop;                % operating (desired) site effect
        a1; a2;             % exponential cbf coefficients
        uMin;               % current minimum input

        % observer (tbd)

        % parameter estimator (tbd)

        % pharmacokinetic parameters (drug)
        k10; k12; k13;      % transport coefficients
        k21; k31; ke0;     
        V1;                 % plasma deposition volume (L)

        % pharmacodynamic parameters (body)
        E0;                 % baseline site effect (BIS)
        Emin;               % desired lower bound of site effect (BIS)
        ce50;               % site concentration at 50% BIS
        gamma;              % cooperativity constant in Hill formula

        % simulation history
        tHist; EHist; xHist; uHist; uMinHist;
        hHist; hDotHist; 

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function sys = PropofolSystem( ...
                k10, k12, k13, k21, k31, ke0, V1, ...
                E0, Emin, ce50, gamma, ...
                a1, a2)

            % initialize system parameters
            sys.k10 = k10; sys.k12 = k12; sys.k13 = k13; 
            sys.k21 = k21; sys.k31 = k31; sys.ke0 = ke0; 
            sys.V1 = V1; sys.E0 = E0; sys.Emin = Emin; 
            sys.ce50 = ce50; sys.gamma = gamma;
            sys.a1 = a1; sys.a2 = a2; sys.uMin = 0; 

            % initialize system matrices and state vectors 
            sys.A = [-k10 - k12 - k13     k12     k13     0
                     k21                   -k21     0       0 
                     k31                    0      -k31     0 
                     ke0                    0       0      -ke0];
            sys.B = [1/V1; 0; 0; 0];
            sys.C = [0 0 0 1];
            sys.E = E0;

            % set lower bound of effect site concentration
            sys.Emin = Emin; 
            sys.ceMin = ce50 * ((E0-Emin)/Emin)^(1/gamma);

            % initialize rows of simulation history vectors
            sys.tHist = zeros(1,1);
            sys.xHist = zeros(4, 1); sys.EHist = zeros(1,1); 
            sys.uHist = zeros(1,1); sys.uMinHist = zeros(1,1);
            sys.hHist = zeros(1, 1); sys.hDotHist = zeros(1, 1);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = computeControl(sys)

            % compute current lower bound of input
            sys.uMin = (sys.V1/sys.ke0)*...
                (sys.ke0 * sys.k10 * sys.x(1) - ...
                sys.ke0 * sys.k12 * (sys.x(2) - sys.x(1)) - ...
                sys.ke0 * sys.k13 * (sys.x(3) - sys.x(1)) - ...
                (sys.a1 + sys.a2) * (sys.ke0 * (sys.x(1) - sys.x(4))) - ...
                sys.a1 * sys.a2 * (sys.x(4) - sys.ceMin));

            % compute control input with switching control policy
            if (sys.uMin > 0)
                sys.u = sys.uMin;
            else
                sys.u = 0;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = updateSiteEffect(sys)

            % compute current site effect
            sys.E = sys.E0 * (1 - ...
                (sys.x(4))^sys.gamma/ ...
                (sys.x(4)^sys.gamma + sys.ce50^sys.gamma));

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = step(sys, dt)

            % compute output ce (in this case a state variable)

            % parameter estimator (tbd)

            % observer (tbd)

            % compute and apply control input
            sys = sys.computeControl(); 

            % plant update
            sys.x = sys.x + dt*(sys.A * sys.x + sys.B * sys.u);

            % site effect update
            sys = sys.updateSiteEffect();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [sys, tHist, xHist, EHist, uHist, uMinHist, ...
                hHist, hDotHist] = simulate(sys, ...
                x0, tStart, tEnd, nTimes)

            % compute discrete time interval
            dt = (tEnd - tStart) / nTimes;

            % reset history vectors
            sys.tHist = zeros(1, nTimes);
            sys.xHist = zeros(4, nTimes); sys.EHist = zeros(1, nTimes); 
            sys.uHist = zeros(1, nTimes); sys.uMinHist = zeros(1, nTimes);
            sys.hHist = zeros(1, nTimes); sys.hDotHist = zeros(1, nTimes);

            % set and log initial conditions
            sys.x = x0; 
            sys.u = 0;
            sys.tHist(1) = tStart;
            sys.xHist(:,1) = x0; 

            for i=2:nTimes

                % time step
                sys = sys.step(dt);

                % log
                sys.tHist(i) = sys.tHist(i-1) + dt;
                sys.xHist(:,i) = sys.x; 
                sys.EHist(i) = sys.E;
                sys.uHist(i) = sys.u;
                sys.uMinHist(i) = sys.uMin;
                sys.hHist(i) = sys.x(4) - sys.ceMin;
                sys.hDotHist(i) = sys.ke0 * (sys.x(1) - sys.x(4));

            end
            
            tHist = sys.tHist; 
            xHist = sys.xHist; EHist = sys.EHist; 
            uHist = sys.uHist; uMinHist = sys.uMinHist;
            hHist = sys.hHist; hDotHist = sys.hDotHist; 

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plotCharacteristicCurve(sys)

            sys.setPlotSettings();     

            % plot characteristic effect curve (not simulated)
            figure()
            ce = linspace(0, 100);
            curve = zeros(100);
            for i=1:100
                curve(i) = sys.E0 * (1 - ...
                ce(i)^sys.gamma /(ce(i)^sys.gamma + sys.ce50^sys.gamma));
            end
            plot(ce, curve);
            xlim([0 6]);
            title('Characteristic propofol effect curve');
            xlabel('Effect site propofol concentration (ug/mL)')
            ylabel('Propofol effect (BIS)')
            ax = gca;
            ax.TitleFontSizeMultiplier = 1.5;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plotInputs(sys)

            sys.setPlotSettings();               

            % plot input and input bound history vs. time
            figure()
            plot(sys.tHist, sys.uHist, 'k'); 
            hold on;
            plot(sys.tHist, sys.uMinHist, 'r');
            legend('Infusion rate', ...
                'Infusion rate lower bound');

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = plotConcentrations(sys)

            sys.setPlotSettings();

            % plot state vector vs. time
            figure()
            plot(sys.tHist, sys.xHist(1,:), 'g');
            hold on
            plot(sys.tHist, sys.xHist(2,:), 'b');
            plot(sys.tHist, sys.xHist(3,:), 'c');
            plot(sys.tHist, sys.xHist(4,:), 'r');
            title('State trajectories in time'); 
            xlabel('Time (min)');
            ylabel('Propofol concentration (ug/mL)')
            legend('Plasma compartment', ...
                'Fast peripheral compartment', ...
                'Slow peripheral compartment', ...
                'Effect site (brain)');
            ax = gca;
            ax.TitleFontSizeMultiplier = 1.5;
            hold off

            % isolate site concentration vs. time
            figure()
            plot(sys.tHist, sys.xHist(4,:), 'r');
            title(['Site propofol concentration (ug/mL)' ...
                ' vs. time (min)']);
            hold on
            yline(sys.ceMin, 'k--', 'lower bound')
            ax = gca;
            ax.TitleFontSizeMultiplier = 1.5;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plotEffect(sys)

            sys.setPlotSettings();

            % plot propofol effect vs. time
            figure()
            plot(sys.tHist, sys.EHist);
            yline(sys.Emin, 'k--', 'lower bound')
            title('Propofol effect (BIS) vs. time (min)');
            ax = gca;
            ax.TitleFontSizeMultiplier = 1.5;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getPhaseData(sys)

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static)

        function setPlotSettings()

            % cosmetic settings
            set(groot, 'defaultTextInterpreter', 'tex');
            set(groot, 'defaultLineLineWidth', 2);
            set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2);  
            set(groot, 'defaultAxesFontSize', 12);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end