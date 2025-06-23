%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define pharmacokinetic coefficients (3-compartment model) 

    % rate of transport into compartment 1 (plasma)
    k10 = 1; k12 = 1; k13 = 1;  
    % infusion deposition volume
    V1 = 1;   
    % rate of transport into compartment 2 (fast peripherals)
    k21 = 1;
    % rate of transport into compartment 3 (slow peripherals)
    k31 = 1; 
    % rate of transport into effect site (brain)
    ke0 = 1; 
    
% define pharmacodynamic coefficients

    % baseline propofol effect with 0 infusion
    E0 = 0; 
    % desired lower bound of propofol effect
    Emin = 20; 
    % site propofol concentration at 50% effect
    ce50 = 20; 
    % cooperativity constant (steepness of Hill function)
    gamma = 2;

% define exponential cbf coefficients
a1 = 1; a2 = 1; 

% define initial conditions
x0 = [0, 0, 0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create system
propofolSystem = propofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                E0, Emin, ce50, gamma, ...
                a1, a2);

% simulate system
tStart = 0; tEnd = 3600; nTimes = 3601;
propofolSystem.simulate(x0, tStart, tEnd, nTimes);

% plot simulation history
propofolSystem.plot();