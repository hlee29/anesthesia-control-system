% runs 12 total simulations of PropofolSystem.m on 1 virtual patient: 
% 6 different initial conditions using exponential controller, 
% 6 different intial conditions using graceful controller. 
% prints safety barriers and state eigenvalues. 
% plots results. 

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
    k10 = 1.79 * 1000 / V1;  % numerator [L min⁻¹] to [mL min⁻¹]
    k12 = 1.75 * 1000 / V2; 
    k13 = 1.11 * 1000 / V3; 
    % transport into compartment 2 (fast peripherals)
    k21 = k12; % symmetric flows for now
    % transport into compartment 3 (slow peripherals)
    k31 = k13; % symmetric flows for now
    % rate of transport into effect site (brain)
    ke0 = 1.24; %[1 min⁻¹]
    
% pharmacodynamic parameters

    % baseline propofol effect with 0 infusion
    E0 = 93;  % [BIS %]
    % desired upper bound of brain activity [BIS %]
    Emax = 60; % 40-60% for surgery
    % site propofol concentration at 50% effect
    ce50 = 3.08; % [ug mL⁻¹]
    % cooperativity constant (steepness of Hill function)
    gamma = 1.68; % source has 1.47 for ce > ce50 and 1.89 for ce < ce50,
                  % took average [dimensionless]

% parameterize observer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

confidence = 0.95;      % confidence level for error bounds (ellipsoid)
P0 = 10 * eye(4);  % initial covariance estimate
Q = 0 * eye(4);  % 4x4 process noise matrix
R = 10^2;                % 1x1 sensor noise matrix (BIS %)
 

% parameterize controller(s) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clamp when uMin > uMax? 
clamp = true; 

% exponential cbf eigenvalues
a1 = 1; a2 = 1; 

% graceful controller coefficients (nonlinear mass-spring-damper)
delta = 1.5; % padding [ug] for first barrier, ceMin + delta
w = 1.9; % natural frequency, must be >0 for negative eigenvalues
z = 1.4; % damping ratio, must be >1 for negative eigenvalues

% actuation limit
uMax = 4 * 10^4; 


% initialize patient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the same male patient with 6 different trajectories (different
% initial conditions)

% for preliminary ceMin computation
x0 = [0; 0; 0; 0]; 
patient0 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x0, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);

x01 = [patient0.ceMin+4; 0; 0; patient0.ceMin+2];
patient1 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x01, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);

x02 = [patient0.ceMin+3; 0; 0; patient0.ceMin+2];
patient2 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x02, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);

x03 = [patient0.ceMin+1; 0; 0; patient0.ceMin+2];
patient3 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x03, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);

x04 = [patient0.ceMin; 0; 0; patient0.ceMin+2.5];
patient4 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x04, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);

x05 = [0; 0; 0; patient0.ceMin+1];
patient5 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x05, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);

x06 = [0; 0; 0; patient0.ceMin+0.3];
patient6 = PropofolSystem(k10, k12, k13, k21, k31, ke0, V1, ...
                x06, E0, Emax, ce50, gamma, ...
                P0, Q, R, confidence, ...
                a1, a2, ...
                uMax, delta, z, w);


% plot these patients' characteristic e-max curve (all same)
patient1.plotCharacteristicCurve();

% print these patients' characteristic eigenvalues (all same)
fprintf('PATIENT INFORMATION: \n')
evals = eig(patient0.A);
fprintf(['Plasma concentration dynamics eigenvalue [min]: ' ...
    '%.4f%+.4fi\n'], real(evals(1)), imag(evals(1)));
fprintf(['Fast peripheral concentration dynamics eigenvalue [min]: ' ...
    '%.4f%+.4fi\n'], real(evals(2)), imag(evals(2)));
fprintf(['Slow peripheral concentration dynamics eigenvalue [min]: ' ...
    '%.4f%+.4fi\n'], real(evals(3)), imag(evals(3)));
fprintf(['Effect site concentration dynamics eigenvalue [min]: ' ...
    '%.4f%+.4fi\n'], real(evals(4)), imag(evals(4)));

% print ceMin corresponding to Emin
fprintf('Absolute concentration barrier ceMin [ug mL⁻¹]: %d\n', ...
    patient1.ceMin);
fprintf('\n');

% print controller clamp setting
if clamp == true fprintf('CLAMPING OF INPUT IS ON.\n');
else fprintf('CLAMPING OF INPUT IS OFF.\n'); end
fprintf('Desired input [ug mL⁻¹ min⁻¹]: %.4\n', patient1.uOp);

fprintf('\n');

% 1. use exponential controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('RUNNING EXPONENTIAL CONTROLLER...\n');

graceful = false;
tStart = 0; tEnd = 60; nTimes = 600;

% scenario 1
patient1.x = x01; xEst0 = [x01(1); x01(2); x01(3); x01(4)];
disp('Running scenario 1...');
tHist = zeros(1, nTimes);
xHist1 = zeros(4, nTimes); EHist1 = zeros(1, nTimes);
uHist1 = zeros(1, nTimes); uMinHist1 = zeros(1, nTimes);
ceDotHist1 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist1, xEstHist1, errorHist1, ...
    EHist1, uHist1, uMinHist1, ...
    ceDotHist1] = patient1.simulate(graceful, clamp, x01, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 2
patient2.x = x02; xEst0 = [x02(1); x02(2); x02(3); x02(4)];
disp('Running scenario 2...');
xHist2 = zeros(4, nTimes); EHist2 = zeros(1, nTimes);
uHist2 = zeros(1, nTimes);uMinHist2 = zeros(1, nTimes);
ceDotHist2 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist2, xEstHist2, errorHist2, ...
    EHist2, uHist2, uMinHist2, ...
    ceDotHist2] = patient2.simulate(graceful, clamp, x02, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 3
patient3.x = x03; xEst0 = [x03(1); x03(2); x03(3); x03(4)];
disp('Running scenario 3...');
xHist3 = zeros(4, nTimes); EHist3 = zeros(1, nTimes);
uHist3 = zeros(1, nTimes); uMinHist3 = zeros(1, nTimes);
ceDotHist3 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist3, xEstHist3, errorHist3, ...
    EHist3, uHist3, uMinHist3, ...
    ceDotHist3] = patient3.simulate(graceful, clamp, x03, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 4
patient4.x = x04; xEst0 = [x04(1); x04(2); x04(3); x04(4)];
disp('Running scenario 4...');
xHist4 = zeros(4, nTimes); EHist4 = zeros(1, nTimes);
uHist4 = zeros(1, nTimes); uMinHist4 = zeros(1, nTimes);
ceDotHist4 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist4, xEstHist4, errorHist4, ...
    EHist4, uHist4, uMinHist4, ...
    ceDotHist4] = patient4.simulate(graceful, clamp, x04, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 5
patient5.x = x05; xEst0 = [x05(1); x05(2); x05(3); x05(4)];
disp('Running scenario 5...');
xHist5 = zeros(4, nTimes); EHist5 = zeros(1, nTimes);
uHist5 = zeros(1, nTimes); uMinHist5 = zeros(1, nTimes);
ceDotHist5 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist5, xEstHist5, errorHist5, ...
    EHist5, uHist5, uMinHist5, ...
    ceDotHist5] = patient5.simulate(graceful, clamp, x05, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 6
patient6.x = x06; xEst0 = [x06(1); x06(2); x06(3); x06(4)];
disp('Running scenario 6...');
xHist6 = zeros(4, nTimes); EHist6 = zeros(1, nTimes);
uHist6 = zeros(1, nTimes); uMinHist6 = zeros(1, nTimes);
ceDotHist6 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist6, xEstHist6, errorHist6, ...
    EHist6, uHist6, uMinHist6, ...
    ceDotHist6] = patient6.simulate(graceful, clamp, x06, xEst0, P0, ...
    tStart, tEnd, nTimes);

fprintf('\n');

% plot effect trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, EHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, EHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, EHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, EHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, EHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, EHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient.Emax, 'k--', 'upper bound for surgery');
yline(40, 'k--', 'lower bound for surgery');
title('Bispectral index trajectories in time (exponential controller)');
xlabel('Time [min]');
ylabel('Bispectral index [%]');

% plot site concentration trajectories
figure('Color', [1 1 1])
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
ylabel('Concentration [ug mL⁻¹]');

% plot inputs
figure('Color', [1 1 1])
hold on
plot(tHist, uHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, uHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, uHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, uHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, uHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, uHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient.uOp, 'k--', 'clinician-set rate');
yline(40000, 'k--', 'maximum rate');
title(['Infusion rate trajectories in time ' ...
    '(exponential controller)']);
xlabel('Time [min]');
ylabel('Infusion rate [ug/(mL * min)]');

% plot phase plane       
figure('Color', [1 1 1])
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
xline(patient.ceMin, 'k--', 'lower bound for surgery', 'LineWidth', 1);  
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0

% plot phase pictures
xlim([-10 10]); ylim([-10 10]);
plot(xHist1(4,:), ceDotHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(x01(4,1), patient.ke0 * (x01(1,1)-x01(4,1)), 'Color', '#66da54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist2(4,:), ceDotHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(x02(4,1), patient.ke0 * (x02(1,1)-x02(4,1)), 'Color','#cada54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist3(4,:), ceDotHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(x03(4,1), patient.ke0 * (x03(1,1)-x03(4,1)), 'Color','#dcd62b', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist4(4,:), ceDotHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(x04(4,1), patient.ke0 * (x04(1,1)-x04(4,1)), 'Color','#dda010', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist5(4,:), ceDotHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(x05(4,1), patient.ke0 * (x05(1,1)-x05(4,1)), 'Color','#dd5510', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist6(4,:), ceDotHist6, 'LineWidth', 2, 'Color', '#980000');
plot(x06(4,1), patient.ke0 * (x06(1,1)-x06(4,1)), 'Color','#980000', ...
    'Marker', '.', 'MarkerSize', 20);
grid 

legend('1st-order safe set', '2nd-order safe set');
title('Safety phase portrait (exponential controller)')
xlabel('Effect site concentration [ug mL⁻¹]')
ylabel('Effect site concentration time derivative [ug mL⁻¹ min⁻¹]')

% 2. simulate graceful controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('RUNNING GRACEFUL CONTROLLER + EKF...\n');

% print ceG, the failsafe barrier
fprintf('Failsafe concentration barrier ceG [ug mL⁻¹]: %d\n', patient.ceG);
graceful = true; % use graceful controller 
tStart = 0; tEnd = 60; nTimes = 600;

% scenario 1
patient1.x = x01; patient1.x = x01; xEst0 = [x01(1); x01(2); x01(3); x01(4)];
patient1.computeErrorBoundEKF(); 
disp('Running scenario 1...');
tHist = zeros(1, nTimes);
xHist1 = zeros(4, nTimes); EHist1 = zeros(1, nTimes);
uHist1 = zeros(1, nTimes); uMinHist1 = zeros(1, nTimes);
ceDotHist1 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist1, xEstHist1, errorHist1, ...
    EHist1, uHist1, uMinHist1, ...
    ceDotHist1] = patient.simulate(graceful, clamp, x01, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 2
patient2.x = x02; xEst0 = [x02(1); x02(2); x02(3); x02(4)];
patient2.computeErrorBoundEKF(); 
disp('Running scenario 2...')
xHist2 = zeros(4, nTimes); EHist2 = zeros(1, nTimes);
uHist2 = zeros(1, nTimes);uMinHist2 = zeros(1, nTimes);
ceDotHist2 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist2, xEstHist2, errorHist2, ...
    EHist2, uHist2, uMinHist2, ...
    ceDotHist2] = patient.simulate(graceful, clamp, x02, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 3
patient3.x = x03; xEst0 = [x03(1); x03(2); x03(3); x03(4)];
patient3.computeErrorBoundEKF();
disp('Running scenario 3...');
xHist3 = zeros(4, nTimes); EHist3 = zeros(1, nTimes);
uHist3 = zeros(1, nTimes); uMinHist3 = zeros(1, nTimes);
ceDotHist3 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist3, xEstHist3, errorHist3, ...
    EHist3, uHist3, uMinHist3, ...
    ceDotHist3] = patient.simulate(graceful, clamp, x03, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 4
patient4.x = x04; xEst0 = [x04(1); x04(2); x04(3); x04(4)];
patient4.computeErrorBoundEKF();
disp('Running scenario 4...');
x04 = [patient.ceMin; 0; 0; patient.ceMin+2.5];
xHist4 = zeros(4, nTimes); EHist4 = zeros(1, nTimes);
uHist4 = zeros(1, nTimes); uMinHist4 = zeros(1, nTimes);
ceDotHist4 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist4, xEstHist4, errorHist4, ...
    EHist4, uHist4, uMinHist4, ...
    ceDotHist4] = patient.simulate(graceful, clamp, x04, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 5
patient5.x = x05; xEst0 = [x05(1); x05(2); x05(3); x05(4)];
patient5.computeErrorBoundEKF();
disp('Running scenario 5...');
x05 = [0; 0; 0; patient.ceMin+1];
xHist5 = zeros(4, nTimes); EHist5 = zeros(1, nTimes);
uHist5 = zeros(1, nTimes); uMinHist5 = zeros(1, nTimes);
ceDotHist5 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist5, xEstHist5, errorHist5, ...
    EHist5, uHist5, uMinHist5, ...
    ceDotHist5] = patient.simulate(graceful, clamp, x05, xEst0, P0, ...
    tStart, tEnd, nTimes);

% scenario 6
patient6.x = x06; 
patient6.computeErrorBoundEKF();
disp('Running scenario 6...');
x06 = [0; 0; 0; patient.ceMin+0.3];
xHist6 = zeros(4, nTimes); EHist6 = zeros(1, nTimes);
uHist6 = zeros(1, nTimes); uMinHist6 = zeros(1, nTimes);
ceDotHist6 = zeros(1, nTimes);
[patient, tHist, tHistContinuous, xHist6, xEstHist6, errorHist6, ...
    EHist6, uHist6, uMinHist6, ...
    ceDotHist6] = patient.simulate(graceful, clamp, x06, xEst0, P0, ...
    tStart, tEnd, nTimes);

% plot effect trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, EHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, EHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, EHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, EHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, EHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, EHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient.Emax, 'k--', 'upper bound for surgery');
yline(40, 'k--', 'lower bound for surgery');
title('Bispectral index trajectories in time (graceful controller)');
xlabel('Time [min]');
ylabel('Bispectral index [%]');

% plot site concentration trajectories
figure('Color', [1 1 1])
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
ylabel('Concentration [ug mL⁻¹]');

% plot inputs
figure('Color', [1 1 1])
hold on
plot(tHist, uHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, uHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, uHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, uHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, uHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, uHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient.uOp, 'k--', 'clinician-set rate');
yline(40000, 'k--', 'maximum rate');
title(['Infusion rate trajectories in time ' ...
    '(graceful controller)']);
xlabel('Time [min]');
ylabel('Infusion rate [ug/(mL * min)]');

% plot phase plane        
figure('Color', [1 1 1])
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
set1 = -1 * lambda1 * (ceVals-patient.ceG)/(patient.ceG-patient.ceMin);
set2 = -1 * lambda2 * (ceVals-patient.ceG)/(patient.ceG-patient.ceMin);
fill([ceVals fliplr(ceVals)], [set1 fliplr(ones(size(ceVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
fill([ceVals fliplr(ceVals)], [set2 fliplr(ones(size(ceVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% plot barriers
plot(ceVals, set1, 'k--', 'LineWidth', 1); 
plot(ceVals, set2, 'k--', 'LineWidth', 1); 
xline(patient.ceMin, 'k--', 'lower bound for surgery', 'LineWidth', 1);
xline(patient.ceG, 'k--', 'failsafe barrier', 'LineWidth', 1); 
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0

% plot 'danger zone' (actuation limits)


% plot phase pictures
xlim([-10 10]); ylim([-10 10]);
plot(xHist1(4,:), ceDotHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(x01(4,1), patient.ke0 * (x01(1,1)-x01(4,1)), 'Color', '#66da54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist2(4,:), ceDotHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(x02(4,1), patient.ke0 * (x02(1,1)-x02(4,1)), 'Color','#cada54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist3(4,:), ceDotHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(x03(4,1), patient.ke0 * (x03(1,1)-x03(4,1)), 'Color','#dcd62b', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist4(4,:), ceDotHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(x04(4,1), patient.ke0 * (x04(1,1)-x04(4,1)), 'Color','#dda010', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist5(4,:), ceDotHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(x05(4,1), patient.ke0 * (x05(1,1)-x05(4,1)), 'Color','#dd5510', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist6(4,:), ceDotHist6, 'LineWidth', 2, 'Color', '#980000');
plot(x06(4,1), patient.ke0 * (x06(1,1)-x06(4,1)), 'Color','#980000', ...
    'Marker', '.', 'MarkerSize', 20);
grid 

legend('Safe set', 'Failsafe set');
title('Safety phase portrait (graceful controller)')
xlabel('Effect site concentration [ug mL⁻¹]')
ylabel('Effect site concentration time derivative [ug mL⁻¹ min⁻¹]')

patient.plotEstimation();