clear all
close all
clc

% parameterize plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters are based on source: eleveld et al., "pkpd model for propofol
% for broad application in anaesthesia and sedation" (2018, british
% journal of anesthesia). units are in [mL] volume, [ug] mass, [min] time. 

% pharmacokinetic parameters

    % deposition volumes
    V1 = 6.28 * 1000; % [L] to [mL]
    V2 = 25.5 * 1000; 
    V3 = 273 * 1000; 
    % transport into/out of compartment 1 (plasma)
    k10 = 1.79 * 1000 / V1;  % numerator [L/min] to [mL/min]
    k12 = 1.75 * 1000 / V2; 
    k13 = 1.11 * 1000 / V3; 
    % transport into compartment 2 (fast peripherals)
    k21 = k12; % symmetric flows for now
    % transport into compartment 3 (slow peripherals)
    k31 = k13; % symmetric flows for now
    % rate of transport into effect site (brain)
    ke0 = 1.24; %[1/min]
    
% pharmacodynamic parameters

    % baseline propofol effect with 0 infusion
    E0 = 93;  % [BIS %]
    % desired lower bound of propofol effect
    Emin = 40; % 40-60% for surgery, 60-80% for sedation
    % site propofol concentration at 50% effect
    ce50 = 3.08; % [ug/mL]
    % cooperativity constant (steepness of Hill function)
    gamma = 1.68; % source has 1.47 for ce > ce50 and 1.89 for ce < ce50,
                  % took average [dimensionless]


% parameterize controller(s) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% smoothening pareto weights
w1 = 0.7; % weight on minimizing input
w2 = 0.3; % weight on minimizing change in input wrt time

% exponential cbf eigenvalues
a1 = 1; a2 = 1; 

% graceful controller coefficients (nonlinear mass-spring-damper)
delta = 1.5; % padding [ug] for first barrier, ceMin + delta
w = 2.1; % natural frequency
z = 1.5; % damping ratio


% 1. simulate exponential controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize patient 
patient = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                E0, Emin, ce50, gamma, ...
                a1, a2, w1, w2, delta, z, w);

% print eigenvalues
evals = eig(patient.A);
fprintf('Plasma concentration dynamics eigenvalue [min]: %.4f%+.4fi\n', real(evals(1)), ...
        imag(evals(1)));
fprintf('Fast peripheral concentration dynamics eigenvalue [min]: %.4f%+.4fi\n', real(evals(2)), ...
        imag(evals(2)));
fprintf('Slow peripheral concentration dynamics eigenvalue [min]: %.4f%+.4fi\n', real(evals(3)), ...
        imag(evals(3)));
fprintf('Effect site concentration dynamics eigenvalue [min]: %.4f%+.4fi\n', real(evals(4)), ...
        imag(evals(4)));

% print ceMin corresponding to Emin
fprintf('Absolute concentration barrier ceMin [ug/mL]: %d\n', patient.ceMin);

% plot patient's characteristic e-max curve
patient.plotCharacteristicCurve();

graceful = false; % use exponential controller 
tStart = 0; tEnd = 60; nTimes = 500;

% scenario 1
x01 = [patient.ceMin+4; 0; 0; patient.ceMin+2];
tHist = zeros(1, nTimes);
xHist1 = zeros(4, nTimes); EHist1 = zeros(1, nTimes);
uHist1 = zeros(1, nTimes); uMinHist1 = zeros(1, nTimes);
ceDotHist1 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist1, EHist1, uHist1, uMinHist1, ...
    ceDotHist1] = patient.simulateContinuous(graceful, x01, tStart, tEnd, nTimes);

% scenario 2
x02 = [patient.ceMin+3; 0; 0; patient.ceMin+2];
xHist2 = zeros(4, nTimes); EHist2 = zeros(1, nTimes);
uHist2 = zeros(1, nTimes);uMinHist2 = zeros(1, nTimes);
ceDotHist2 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist2, EHist2, uHist2, uMinHist2, ...
    ceDotHist2] = patient.simulateContinuous(graceful, x02, tStart, tEnd, nTimes);

% scenario 3
x03 = [patient.ceMin+1; 0; 0; patient.ceMin+2];
xHist3 = zeros(4, nTimes); EHist3 = zeros(1, nTimes);
uHist3 = zeros(1, nTimes); uMinHist3 = zeros(1, nTimes);
ceDotHist3 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist3, EHist3, uHist3, uMinHist3, ...
    ceDotHist3] = patient.simulateContinuous(graceful, x03, tStart, tEnd, nTimes);

% scenario 4
x04 = [patient.ceMin; 0; 0; patient.ceMin+2.5];
xHist4 = zeros(4, nTimes); EHist4 = zeros(1, nTimes);
uHist4 = zeros(1, nTimes); uMinHist4 = zeros(1, nTimes);
ceDotHist4 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist4, EHist4, uHist4, uMinHist4, ...
    ceDotHist4] = patient.simulateContinuous(graceful, x04, tStart, tEnd, nTimes);

% scenario 5
x05 = [0; 0; 0; patient.ceMin+1];
xHist5 = zeros(4, nTimes); EHist5 = zeros(1, nTimes);
uHist5 = zeros(1, nTimes); uMinHist5 = zeros(1, nTimes);
ceDotHist5 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist5, EHist5, uHist5, uMinHist5, ...
    ceDotHist5] = patient.simulateContinuous(graceful, x05, tStart, tEnd, nTimes);

% scenario 6
x06 = [0; 0; 0; patient.ceMin+0.3];
xHist6 = zeros(4, nTimes); EHist6 = zeros(1, nTimes);
uHist6 = zeros(1, nTimes); uMinHist6 = zeros(1, nTimes);
ceDotHist6 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist6, EHist6, uHist6, uMinHist6, ...
    ceDotHist6] = patient.simulateContinuous(graceful, x06, tStart, tEnd, nTimes);

% plot effect trajectories
figure()
hold on
plot(tHist, EHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, EHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, EHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, EHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, EHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, EHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient.Emin, 'k--', 'lower bound for surgery');
title('Propofol effect trajectories in time (exponential controller)');
xlabel('Time [min]');
ylabel('Total propofol effect [decrease in BIS %]');

% plot site concentration trajectories
figure()
hold on
plot(tHist, xHist1(4,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xHist2(4,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xHist3(4,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xHist4(4,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xHist5(4,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xHist6(4,:), 'LineWidth', 2, 'Color', '#980000');
yline(patient.ceMin, 'k--', 'lower bound for surgery');
title(['Effect site concentration trajectories in time ' ...
    '(exponential controller)']);
xlabel('Time [min]');
ylabel('Concentration [ug/mL]');

% plot inputs
figure()
hold on
plot(tHist, uHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, uHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, uHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, uHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, uHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, uHist6, 'LineWidth', 2, 'Color', '#980000');
title(['Infusion rate trajectories in time ' ...
    '(exponential controller)']);
xlabel('Time [min]');
ylabel('Infusion rate [ug/(mL * min)]');

% plot phase plane
        
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

% plot barriers
plot(ceVals, bound, 'k--', 'LineWidth', 1); 
xline(patient.ceMin, 'k--', 'lower bound for surgery', 'LineWidth', 1);  % h = 0
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0

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
plot(xHist6(4,:), ceDotHist6, 'LineWidth', 2, 'Color', '#980000');
plot(x06(4,1), patient.ke0 * (x06(1,1)-x06(4,1)), 'Color','#980000', 'Marker', '.', 'MarkerSize', 20);
grid 

legend('1st-order safe set', '2nd-order safe set');
title('Safety phase portrait (exponential controller)')
xlabel('Effect site concentration [ug/mL]')
ylabel('Effect site concentration time derivative [ug/(mL*s)]')

% 2. simulate graceful controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print ceG, the failsafe barrier
fprintf('Failsafe concentration barrier ceG [ug/mL]: %d\n', patient.ceG);
graceful = true; % use graceful controller 
tStart = 0; tEnd = 60; nTimes = 800;

% scenario 1
x01 = [patient.ceMin+4; 0; 0; patient.ceMin+2];
tHist = zeros(1, nTimes);
xHist1 = zeros(4, nTimes); EHist1 = zeros(1, nTimes);
uHist1 = zeros(1, nTimes); uMinHist1 = zeros(1, nTimes);
ceDotHist1 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist1, EHist1, uHist1, uMinHist1, ...
    ceDotHist1] = patient.simulateContinuous(graceful, x01, tStart, tEnd, nTimes);

% scenario 2
x02 = [patient.ceMin+3; 0; 0; patient.ceMin+2];
xHist2 = zeros(4, nTimes); EHist2 = zeros(1, nTimes);
uHist2 = zeros(1, nTimes);uMinHist2 = zeros(1, nTimes);
ceDotHist2 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist2, EHist2, uHist2, uMinHist2, ...
    ceDotHist2] = patient.simulateContinuous(graceful, x02, tStart, tEnd, nTimes);

% scenario 3
x03 = [patient.ceMin+1; 0; 0; patient.ceMin+2];
xHist3 = zeros(4, nTimes); EHist3 = zeros(1, nTimes);
uHist3 = zeros(1, nTimes); uMinHist3 = zeros(1, nTimes);
ceDotHist3 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist3, EHist3, uHist3, uMinHist3, ...
    ceDotHist3] = patient.simulateContinuous(graceful, x03, tStart, tEnd, nTimes);

% scenario 4
x04 = [patient.ceMin; 0; 0; patient.ceMin+2.5];
xHist4 = zeros(4, nTimes); EHist4 = zeros(1, nTimes);
uHist4 = zeros(1, nTimes); uMinHist4 = zeros(1, nTimes);
ceDotHist4 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist4, EHist4, uHist4, uMinHist4, ...
    ceDotHist4] = patient.simulateContinuous(graceful, x04, tStart, tEnd, nTimes);

% scenario 5
x05 = [0; 0; 0; patient.ceMin+1];
xHist5 = zeros(4, nTimes); EHist5 = zeros(1, nTimes);
uHist5 = zeros(1, nTimes); uMinHist5 = zeros(1, nTimes);
ceDotHist5 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist5, EHist5, uHist5, uMinHist5, ...
    ceDotHist5] = patient.simulateContinuous(graceful, x05, tStart, tEnd, nTimes);

% scenario 6
x06 = [0; 0; 0; patient.ceMin+0.3];
xHist6 = zeros(4, nTimes); EHist6 = zeros(1, nTimes);
uHist6 = zeros(1, nTimes); uMinHist6 = zeros(1, nTimes);
ceDotHist6 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist6, EHist6, uHist6, uMinHist6, ...
    ceDotHist6] = patient.simulateContinuous(graceful, x06, tStart, tEnd, nTimes);

% plot effect trajectories
figure()
hold on
plot(tHist, EHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, EHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, EHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, EHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, EHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, EHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient.Emin, 'k--', 'lower bound for surgery');
title('Propofol effect trajectories in time (graceful controller)');
xlabel('Time [min]');
ylabel('Total propofol effect [decrease in BIS %]');

% plot site concentration trajectories
figure()
hold on
plot(tHist, xHist1(4,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xHist2(4,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xHist3(4,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xHist4(4,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xHist5(4,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xHist6(4,:), 'LineWidth', 2, 'Color', '#980000');
yline(patient.ceG, 'k--', 'failsafe barrier');
yline(patient.ceMin, 'k--', 'lower bound for surgery');
title(['Effect site concentration trajectories in time ' ...
    '(graceful controller)']);
xlabel('Time [min]');
ylabel('Concentration [ug/mL]');

% plot inputs
figure()
hold on
plot(tHist, uHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, uHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, uHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, uHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, uHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, uHist6, 'LineWidth', 2, 'Color', '#980000');
title(['Infusion rate trajectories in time ' ...
    '(graceful controller)']);
xlabel('Time [min]');
ylabel('Infusion rate [ug/(mL * min)]');

% plot phase plane        
figure()
hold on

% shade h >= 0 
fill([patient.ceG 100 100 patient.ceG], [-50 -50 50 50], [0.4 0.4 0.4], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.6);
fill([patient.ceMin 100 100 patient.ceMin], [-50 -50 50 50], [0.7 0.7 0.7], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);

% shade new secondary barrier
ceVals = linspace(-100, 100, 500);
lambda1 = -1 * patient.w * sqrt(patient.z^2 - 1) + patient.z*patient.w; 
lambda2 = patient.w * sqrt(patient.z^2 - 1) + patient.z*patient.w; 
bound1 = -1 * lambda1 * (ceVals-patient.ceG)/(patient.ceG-patient.ceMin);
bound2 = -1 * lambda2 * (ceVals-patient.ceG)/(patient.ceG-patient.ceMin);
fill([ceVals fliplr(ceVals)], [bound1 fliplr(ones(size(ceVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
fill([ceVals fliplr(ceVals)], [bound2 fliplr(ones(size(ceVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% plot barriers
plot(ceVals, bound1, 'k--', 'LineWidth', 1, 'DisplayName', 'eigenvalue set'); 
plot(ceVals, bound2, 'k--', 'LineWidth', 1, 'DisplayName', 'eigenvalue set'); 
xline(patient.ceMin, 'k--', 'lower bound for surgery', 'LineWidth', 1);  % h = 0
xline(patient.ceG, 'k--', 'failsafe barrier', 'LineWidth', 1); 
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0


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
plot(xHist6(4,:), ceDotHist6, 'LineWidth', 2, 'Color', '#980000');
plot(x06(4,1), patient.ke0 * (x06(1,1)-x06(4,1)), 'Color','#980000', 'Marker', '.', 'MarkerSize', 20);
grid 

legend('Safe set', 'Failsafe set');
title('Safety phase portrait (graceful controller)')
xlabel('Effect site concentration [ug/mL]')
ylabel('Effect site concentration time derivative [ug/(mL*s)]')