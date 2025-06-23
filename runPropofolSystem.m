clear all
close all
clc

% parameterize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define pharmacokinetic coefficients (3-compartment model) 
% numbers from ref. 38. use [mL], [ug], [min]. 

    % deposition volumes
    V1 = 6.28 * 1000; % (L > mL)
    V2 = 25.5 * 1000; 
    V3 = 273 * 1000; 
    % transport into/out of compartment 1 (plasma)
    k10 = 1.79 * 1000 / V1;  % numerator L/min > mL/min
    k12 = 1.75 * 1000 / V2; 
    k13 = 1.11 * 1000 / V3; 
    % transport into compartment 2 (fast peripherals)
    k21 = k12; % symmetric flows for now
    % transport into compartment 3 (slow peripherals)
    k31 = k13; % symmetric flows for now
    % rate of transport into effect site (brain)
    ke0 = 1.24; %(1/min) 
    
% define pharmacodynamic coefficients

    % baseline propofol effect with 0 infusion
    E0 = 93;  % (BIS %)
    % desired lower bound of propofol effect
    Emin = 40; % google: 40-60% for surgery, 60-80% for sedation
    % site propofol concentration at 50% effect
    ce50 = 3.08; % (ug/mL)
    % cooperativity constant (steepness of Hill function)
    gamma = 1.68; % ref. 38 has 1.47 for ce > ce50 and 1.89 for ce < ce50,
                  % took average. (dimensionless)

% define exponential cbf coefficients
a1 = 1; a2 = 1; 

% 1. start in safe set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
patient = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                E0, Emin, ce50, gamma, ...
                a1, a2);

% plot characteristic E-max curve
patient.plotCharacteristicCurve();

tStart = 0; tEnd = 30; nTimes = 200;

% scenario 1
x0 = [patient.ceMin+4; 0; 0; patient.ceMin+2];
tStart = 0; tEnd = 30; nTimes = 100;
tHist = zeros(1, nTimes);
xHist1 = zeros(4, nTimes); EHist1 = zeros(1, nTimes);
uHist1 = zeros(1, nTimes); uMinHist1 = zeros(1, nTimes);
hHist1 = zeros(1, nTimes); hDotHist1 = zeros(1, nTimes);
[patient, tHist, xHist1, EHist1, uHist1, uMinHist1, ...
    hHist1, hDotHist1] = patient.simulate(x0, tStart, tEnd, nTimes);

% scenario 2
x0 = [patient.ceMin+2; 0; 0; patient.ceMin+2];
xHist2 = zeros(4, nTimes); EHist2 = zeros(1, nTimes);
uHist2 = zeros(1, nTimes);uMinHist2 = zeros(1, nTimes);
hHist2 = zeros(1, nTimes); hDotHist2 = zeros(1, nTimes);
[patient, tHist, xHist2, EHist2, uHist2, uMinHist2, ...
    hHist2, hDotHist2] = patient.simulate(x0, tStart, tEnd, nTimes);

% scenario 3
x0 = [patient.ceMin+1; 0; 0; patient.ceMin+2];
xHist3 = zeros(4, nTimes); EHist3 = zeros(1, nTimes);
uHist3 = zeros(1, nTimes); uMinHist3 = zeros(1, nTimes);
hHist3 = zeros(1, nTimes); hDotHist3 = zeros(1, nTimes);
[patient, tHist, xHist3, EHist3, uHist3, uMinHist3, ...
    hHist3, hDotHist3] = patient.simulate(x0, tStart, tEnd, nTimes);

% scenario 4
x0 = [patient.ceMin; 0; 0; patient.ceMin+2];
xHist4 = zeros(4, nTimes); EHist4 = zeros(1, nTimes);
uHist4 = zeros(1, nTimes); uMinHist4 = zeros(1, nTimes);
hHist4 = zeros(1, nTimes); hDotHist4 = zeros(1, nTimes);
[patient, tHist, xHist4, EHist4, uHist4, uMinHist4, ...
    hHist4, hDotHist4] = patient.simulate(x0, tStart, tEnd, nTimes);

% scenario 5
x0 = [patient.ceMin-4; 0; 0; patient.ceMin+4];
xHist5 = zeros(4, nTimes); EHist5 = zeros(1, nTimes);
uHist5 = zeros(1, nTimes); uMinHist5 = zeros(1, nTimes);
hHist5 = zeros(1, nTimes); hDotHist5 = zeros(1, nTimes);
[patient, tHist, xHist5, EHist5, uHist5, uMinHist5, ...
    hHist5, hDotHist5] = patient.simulate(x0, tStart, tEnd, nTimes);

% plots and more plots%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
figure()
hold on
% shade h >= 0 
fill([0 10 10 0], [-10 -10 10 10], [0.9 0.9 0.9], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.3);
            
% shade hDot >= -a1 * h
hVals = linspace(-10, 10, 500);
fill([hVals fliplr(hVals)], [-hVals fliplr(ones(size(hVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% plot phase pictures
xlim([-10 10]); ylim([-10 10]);
plot(hHist1, hDotHist1, 'LineWidth', 1.5, 'Color', '#66da54');
plot(hHist2, hDotHist2, 'LineWidth', 1.5, 'Color', '#cada54');
plot(hHist3, hDotHist3, 'LineWidth', 1.5, 'Color', '#dcd62b');
plot(hHist4, hDotHist4, 'LineWidth', 1.5, 'Color', '#dda010');
plot(hHist5, hDotHist5, 'LineWidth', 1.5, 'Color', '#dd5510');

% plot barriers
plot(hVals, -patient.a1 * hVals, 'k--', 'LineWidth', 1.5); 
xline(0, 'k--');  % h = 0

legend('dh(0)/dt >= 0', 'dh(0)/dt < 0');
title('Safety phase plane')
xlabel('Barrier function h')
ylabel('Barrier function time derivative dh/dt')
