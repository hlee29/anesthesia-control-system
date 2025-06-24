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

% define minimization smoothening weights
w1 = 1; w2 = 0;

% 1. start in safe set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
patient = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                E0, Emin, ce50, gamma, ...
                a1, a2, w1, w2);
fprintf('ceMin = %d\n', patient.ceMin);

% plot characteristic E-max curve
patient.plotCharacteristicCurve();

tStart = 0; tEnd = 60; nTimes = 500;

% scenario 1
x01 = [patient.ceMin+4; 0; 0; patient.ceMin+2];
tStart = 0; tEnd = 30; nTimes = 100;
tHist = zeros(1, nTimes);
xHist1 = zeros(4, nTimes); EHist1 = zeros(1, nTimes);
uHist1 = zeros(1, nTimes); uMinHist1 = zeros(1, nTimes);
ceDotHist1 = zeros(1, nTimes);
[patient, tHist, xHist1, EHist1, uHist1, uMinHist1, ...
    ceDotHist1] = patient.simulate(x01, tStart, tEnd, nTimes);

% scenario 2
x02 = [patient.ceMin+2; 0; 0; patient.ceMin+2];
xHist2 = zeros(4, nTimes); EHist2 = zeros(1, nTimes);
uHist2 = zeros(1, nTimes);uMinHist2 = zeros(1, nTimes);
ceDotHist2 = zeros(1, nTimes);
[patient, tHist, xHist2, EHist2, uHist2, uMinHist2, ...
    ceDotHist2] = patient.simulate(x02, tStart, tEnd, nTimes);

% scenario 3
x03 = [patient.ceMin+1; 0; 0; patient.ceMin+2];
xHist3 = zeros(4, nTimes); EHist3 = zeros(1, nTimes);
uHist3 = zeros(1, nTimes); uMinHist3 = zeros(1, nTimes);
ceDotHist3 = zeros(1, nTimes);
[patient, tHist, xHist3, EHist3, uHist3, uMinHist3, ...
    ceDotHist3] = patient.simulate(x03, tStart, tEnd, nTimes);

% scenario 4
x04 = [patient.ceMin; 0; 0; patient.ceMin+2];
xHist4 = zeros(4, nTimes); EHist4 = zeros(1, nTimes);
uHist4 = zeros(1, nTimes); uMinHist4 = zeros(1, nTimes);
ceDotHist4 = zeros(1, nTimes);
[patient, tHist, xHist4, EHist4, uHist4, uMinHist4, ...
    ceDotHist4] = patient.simulate(x04, tStart, tEnd, nTimes);

% scenario 5
x05 = [0; 0; 0; patient.ceMin+1];
xHist5 = zeros(4, nTimes); EHist5 = zeros(1, nTimes);
uHist5 = zeros(1, nTimes); uMinHist5 = zeros(1, nTimes);
ceDotHist5 = zeros(1, nTimes);
[patient, tHist, xHist5, EHist5, uHist5, uMinHist5, ...
    ceDotHist5] = patient.simulate(x05, tStart, tEnd, nTimes);

patient.plotInputs();

% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
hold on
plot(tHist, EHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, EHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, EHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, EHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, EHist5, 'LineWidth', 2, 'Color', '#dd5510');
yline(patient.Emin, 'k--', 'lower bound');
title('Propofol effect trajectories in time');
xlabel('Time [min]');
ylabel('Total propofol effect [decrease in BIS %]');

figure()
hold on
plot(tHist, xHist1(4,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xHist2(4,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xHist3(4,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xHist4(4,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xHist5(4,:), 'LineWidth', 2, 'Color', '#dd5510');
yline(patient.ceMin, 'k--', 'lower bound');
title('Effect site concentration trajectories in time');
xlabel('Time [min]');
ylabel('Concentration [ug/mL]');

% plot phase plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
figure()
hold on
% shade h >= 0 
fill([patient.ceMin 20 20 patient.ceMin], [-10 -10 10 10], [0.9 0.9 0.9], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);
            
% shade hDot >= -a1 * h
ceVals = linspace(-20, 20, 500);
bound = -1 * patient.a1 * (ceVals - patient.ceMin);
fill([ceVals fliplr(ceVals)], [bound fliplr(ones(size(ceVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% plot phase pictures
xlim([-10 10]); ylim([-10 10]);
plot(xHist1(4,:), ceDotHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(x01(4,1), patient.ke0 * (x01(1,1)-x01(4,1)), 'Color', '#66da54', 'Marker', '.', 'MarkerSize', 20);
plot(xHist2(4,:), ceDotHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(x02(4,1), patient.ke0 * (x02(1,1)-x02(4,1)), 'Color','#cada54', 'Marker', '.', 'MarkerSize', 20);
plot(xHist3(4,:), ceDotHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(x03(4,1), patient.ke0 * (x03(1,1)-x03(4,1)), 'Color','#dcd62b', 'Marker', '.', 'MarkerSize', 20);
plot(xHist4(4,:), ceDotHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(x04(4,1), patient.ke0 * (x04(1,1)-x04(4,1)), 'Color','#dda010', 'Marker', '.', 'MarkerSize', 20);
plot(xHist5(4,:), ceDotHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(x05(4,1), patient.ke0 * (x05(1,1)-x05(4,1)), 'Color','#dd5510', 'Marker', '.', 'MarkerSize', 20);
grid

% plot barriers
plot(ceVals, bound, 'k--', 'LineWidth', 1); 
xline(patient.ceMin, 'k--', 'LineWidth', 1);  % h = 0
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0


legend('1st-order safe set', '2nd-order safe set');
title('Safety phase plane')
xlabel('Effect site concentration [ug/mL]')
ylabel('Effect site concentration time derivative [ug/(mL*s)]')
