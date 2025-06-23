classdef propofolSystem

    % an automated IV propofol infusion control system

    properties

        % linear plant: infusion mass/vol input, site concentration output
        A; B; C;            % system matrices
        x;                  % current system state
        u;                  % current input (mg/mL*s)
        E;                  % current site effect (BIS)
        ceDot;              % time derivative of ce

        % optimal controller
        ceMin;              % site concentration lower bound
        Eop;                % operating (desired) site effect
        a1; a2;             % exponential cbf coefficients

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
        tHist; EHist; xHist; uHist;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function sys = propofolSystem( ...
                k10, k12, k13, k21, k31, ke0, V1, ...
                E0, Emin, ce50, gamma, ...
                a1, a2)

            % initialize system parameters
            sys.k10 = k10; sys.k12 = k12; sys.k13 = k13; 
            sys.k21 = k21; sys.k31 = k31; sys.ke0 = ke0; 
            sys.V1 = V1; sys.E0 = E0; sys.Emin = Emin; 
            sys.ce50 = ce50; sys.gamma = gamma;
            sys.a1 = a1; sys.a2 = a2; 

            % initialize system matrices and state vectors 
            sys.A = [-(k10 + k12 + k13), k12, k13, 0;
                     k21, -k21, 0, 0; 
                     k31, 0, -k31, 0; 
                     ke0, 0, 0, -ke0];
            sys.B = [1/V1, 0, 0, 0];
            sys.C = [1, 0, 0, 0];
            sys.E = E0;

            % set lower bound of effect site concentration
            sys.Emin = Emin; 
            sys.ceMin = ce50 * ((E0-Emin)/Emin)^(1/gamma);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = computeControl(sys)

            % compute current lower bound of input
            uMin = (sys.V1/sys.ke0)*...
                (sys.ke0 * sys.k10 * sys.x(1) - ...
                sys.ke0 * sys.k12 * (sys.x(2) - sys.x(1)) - ...
                sys.ke0 * sys.k13 * (sys.x(3) - sys.x(1)) - ...
                (sys.a1 + sys.a2) * sys.ceDot - ...
                sys.a1 * sys.a2 * (sys.x(4) - sys.ceMin));

            % compute control input with switching control policy
            if (uMin > 0)
                sys.u = uMin;
            else
                sys.u = 0;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = updateSiteEffect(sys)

            % compute current site effect
            sys.E = sys.E0 * (1 - ...
                (sys.ceMin)^sys.gamma/ ...
                (sys.ceMin^sys.gamma + sys.ce50^sys.gamma));

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = step(sys, dt)

            % compute output ce (in this case a state variable)

            % parameter estimator (tbd)

            % observer (tbd)

            % compute control
            sys.u = computeControl(); 

            % plant update
            sys.x = sys.x + dt*(sys.A * sys.x + sys.B * sys.u);

            % site effect update
            sys.updateSiteEffect();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = simulate(sys, ...
                x0, tStart, tEnd, nTimes)

            % compute discrete time interval
            dt = (tStart - tEnd) / nTimes;

            % reset history vectors
            sys.xHist = zeros(4, nTimes);
            sys.EHist = zeroes(1, nTimes); 
            sys.uHist = zeroes(1, nTimes); 
            sys.tHist = zeroes(1, nTimes);

            % set and log initial conditions
            sys.x = x0; 
            sys.xHist(:,1) = x0; 

            for i=2:nTimes

                % time step
                sys = step(sys, dt);

                % log
                sys.tHist(i) = sys.tHist(i-1) + dt;
                sys.xHist(:,i) = sys.x; 
                sys.EHist(i) = sys.E;
                sys.uHist(i) = sys.u;

            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sys = plot(sys)

            close all
            set(0,'defaulttextinterpreter','latex')
            set(groot,'defaultLineLineWidth',1.5)
            set(0, 'DefaultFigureColor', [1 1 1])

            % plot propofol effect vs. time
            figure()
            plot(sys.tHist, sys.EHist);
            title('Propofol effect (BIS) vs. time (s)')

            % plot state vector vs. time
            figure()
            plot(sys.tHist, sys.xHist(1,:), 'g', ...
                sys.tHist, sys.xHist(2,:), 'b', ...
                sys.tHist, sys.xHist(3,:), 'c', ...
                sys.tHist, sys.xHist(4,:), 'r');
            title('Compartment propofol concentrations (mg/mL vs. time (s)');
            legend('Plasma concentration', ...
                'Fast peripheral concentration', ...
                'Slow peripheral concentration', ...
                'Site (brain) concentration');

            % isolate site concentration vs. time
            figure()
            plot(sys.tHist, sys.xHist(4,:), 'r');
            title('Site propofol concentration (mg/mL vs. time (s)')

            % plot characteristic effect curve (not simulated)
            figure()
            title(['Characteristic curve: ' ...
                'Propofol effect (BIS) vs. ' ...
                'Site propofol concentration (mg/mL'])
            
            curve = @(x) sys.E0 * (1 - ...
                x^sys.gamma /(x^sys.gamma + sys.ce50^sys.gamma));
            fplot(curve, [0, 5]);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end