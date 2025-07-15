% runs 12 total simulations of PropofolSystem.m on 1 virtual patient: 
% 6 different initial conditions using exponential controller, 
% 6 different intial conditions using graceful controller. 
% prints safety barriers and state eigenvalues. 
% plots results. 

clear all
close all
clc

set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

% parameterize plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters are based on source: eleveld et al., "pkpd model for
% propofol for broad application in anaesthesia and sedation" (2018,
% british journal of anesthesia).
% units are in [mL] volume, [ug] mass, [min] time. 

% pharmacokinetic parameters

    % deposition volumes
    V1 = 6.28 * 1000; % [L] to [mL]
    V2 = 25.5 * 1000; 
    V3 = 273 * 1000; 
    % transport into/out of compartment 1 (plasma)
    k10 = 1.79 * 1000 / V1;  % numerator [L min$^{-1}$] to [mL min$^{-1}$]
    k12 = 1.75 * 1000 / V2; 
    k13 = 1.11 * 1000 / V3; 
    % transport into compartment 2 (fast peripherals)
    k21 = k12; % symmetric flows for now
    % transport into compartment 3 (slow peripherals)
    k31 = k13; % symmetric flows for now
    % rate of transport into effect site (brain)
    ke0 = 1.24; %[1 min$^{-1}$]
    
% pharmacodynamic parameters

    % baseline propofol effect with 0 infusion
    BIS0 = 93;  % [BIS %]
    % site propofol concentration at 50% effect
    ce50 = 3.08; % [$\mu$g mL$^{-1}$]
    % cooperativity constant (steepness of Hill function)
    gamma = 1.68; % source has 1.47 for ce > ce50 and 1.89 for ce < ce50,
                  % took average [dimensionless]

% parameterize observer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

observer = false;        % use observer? 
confidence = 0.95;       % confidence level for error bounds (ellipsoid)
P0 = 10 * eye(4);        % initial covariance estimate
Q = 10^(-2) * eye(4);    % 4x4 process noise matrix
R = 100;                 % 1x1 sensor noise matrix
 

% parameterize controllers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general settings
uMax = 5 * 10^4;        % maximum bolus dosage
BISmax = 60;            % maximum for surgery 
BISdes = 50;            % good for surgery
clamp = true;          % clamp when uMin > uMax?

% exponential cbf eigenvalues
a1 = 1; a2 = 1; 

% graceful controller coefficients (nonlinear mass-spring-damper)
delta = 1.5;            % padding [ug] for first barrier, ceMin + delta
w = 1.9;                % natural freq., must be > 0 for eigenvalues < 0
z = 1.2;                % damping ratio, must be >1 for eigenvalues < 0 


% initialize patient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a 70 kg adult male patient undergoing 6 different trajectories

% for preliminary ceMin computation
x0 = [0; 0; 0; 0]; xEst0 = x0; 
patient0 = PropofolSystem(x0, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst0, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);
patient0.plotCharacteristicCurve(); % plot patient's Emax curve

x01 = [patient0.ceMin+4; 0; 0; patient0.ceMin+2]; xEst01 = x01;
patient1 = PropofolSystem(x01, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst01, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);

x02 = [patient0.ceMin+3; 0; 0; patient0.ceMin+2]; xEst02 = x02; 
patient2 = PropofolSystem(x02, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst02, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);

x03 = [patient0.ceMin+1; 0; 0; patient0.ceMin+2]; xEst03 = x03; 
patient3 = PropofolSystem(x03, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst03, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);

x04 = [patient0.ceMin; 0; 0; patient0.ceMin+2.5]; xEst04 = x04; 
patient4 = PropofolSystem(x04, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst04, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);

x05 = [0; 0; 0; patient0.ceMin+1]; xEst05 = x05; 
patient5 = PropofolSystem(x05, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst05, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);

x06 = [0; 0; 0; patient0.ceMin+0.3]; xEst06 = x06; 
patient6 = PropofolSystem(x06, k10, k12, k13, k21, k31, ke0, V1, ...
                BIS0, BISdes, BISmax, ce50, gamma, ...
                xEst06, P0, Q, R, ...
                uMax, a1, a2, delta, z, w);

% print this patient's characteristic eigenvalues (all same)
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

% print ceMin corresponding to BISmax
fprintf("Absolute concentration barrier ceMin [$\mu$g mL$^{-1}$]: %d\n", ...
    patient1.ceMin);
fprintf('\n');

% print controller settings
fprintf('CONTROLLER SETTINGS:\n')
fprintf('Exponential controller: a1 = %.2f, a2 = %.2f', a1, a2);
fprintf(['Graceful controller: omega = %.2f, zeta = %.2f, ' ...
    'failsafe barrier ceG [$\mu$g mL$^{-1}$] = %.4f\n'], w, z, patient0.ceG);
if clamp == true fprintf('Clamping of input is ON.\n');
else fprintf('Clamping of input is OFF.\n'); end
fprintf('\n');

% print observer settings
H = [0 0 0 -1 * patient0.BIS0 * patient0.gamma * ...
                (patient0.ce50^(patient0.gamma) * ...
                patient0.ceDes^(patient0.gamma-1)) / ...
                (patient0.ceDes^(patient0.gamma) + ...
                patient0.ce50^(patient0.gamma))^2];
obs = [H ; H * patient0.A; H * patient0.A^2; H * patient0.A^3]; 
obsRank = rank(obs);
singVals = svd(obs);
conditioningNum = max(singVals) / min(singVals); 
fprintf('OBSERVER SETTINGS:\n')
fprintf(['Rank of observability matrix around desired ' ...
    'concentration = %d\n'], obsRank);
fprintf('Singular value of c1 = %d, ', singVals(1));
fprintf('singular value of c2 = %d;\n', singVals(2));
fprintf('Singular value of c3 = %d, ', singVals(3));
fprintf('singular value of ce = %d\n', singVals(4));
fprintf('Conditioning number of observability matrix = %d\n', ...
    conditioningNum);
if observer == true 
    fprintf('Controller and estimator are CONNECTED.\n')
else 
    fprintf('Controller and estimator are DISCONNECTED.\n')
end
fprintf('\n');


% 1. use exponential controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('RUNNING EXPONENTIAL CONTROLLER...\n');

graceful = false;
tStart = 0; tEnd = 60; nTimes = 600;
tHist = zeros(1, nTimes);

% scenario 1
disp('Running scenario 1...');
[patient1, tHist, tHistContinuous1, xHist1, xEstHist1,...
    uHist1, uMinHist1, ceDotHist1,...
    BISHist1, BISEstHist1] = patient1.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 2
disp('Running scenario 2...');
[patient2, tHist2, tHistContinuous2, xHist2, xEstHist2,...
    uHist2, uMinHist2, ceDotHist2,...
    BISHist2, BISEstHist2] = patient2.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 3
disp('Running scenario 3...');
[patient3, tHist3, tHistContinuous3, xHist3, xEstHist3,...
    uHist3, uMinHist3, ceDotHist3,...
    BISHist3, BISEstHist3] = patient3.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 4
disp('Running scenario 4...');
[patient4, tHist4, tHistContinuous4, xHist4, xEstHist4,...
    uHist4, uMinHist4, ceDotHist4,...
    BISHist4, BISEstHist4] = patient4.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 5
disp('Running scenario 5...');
[patient5, tHist5, tHistContinuous5, xHist5, xEstHist5,...
    uHist5, uMinHist5, ceDotHist5,...
    BISHist5, BISEstHist5] = patient5.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 6
disp('Running scenario 6...');
[patient6, tHist6, tHistContinuous6, xHist6, xEstHist6,...
    uHist6, uMinHist6, ceDotHist6,...
    BISHist6, BISEstHist6] = patient6.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

fprintf('\n');

% plot effect trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, BISHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, BISHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, BISHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, BISHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, BISHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, BISHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient0.BISmax, 'k--', 'upper bound for surgery');
yline(patient0.BISdes, 'k--', 'desired');
yline(40, 'k--', 'lower bound for surgery');
title('Bispectral index trajectories in time (exponential controller)');
xlabel('Time [min]');
ylabel('Bispectral index [\%]');

% plot site concentration trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, xHist1(4,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xHist2(4,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xHist3(4,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xHist4(4,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xHist5(4,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xHist6(4,:), 'LineWidth', 2, 'Color', '#980000');
yline(patient0.ceMin, 'k--', 'lower bound for surgery');
title(['Effect site concentration trajectories in time ' ...
    '(exponential controller)']);
xlabel('Time [min]');
ylabel('Concentration [$\mu$g mL$^{-1}$]');

% plot inputs
figure('Color', [1 1 1])
hold on
plot(tHist, uHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, uHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, uHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, uHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, uHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, uHist6, 'LineWidth', 2, 'Color', '#980000');
yline(uMax, 'k--', 'maximum rate');
title(['Infusion rate trajectories in time ' ...
    '(exponential controller)']);
xlabel('Time [min]');
ylabel('Infusion rate [$\mu$g mL$^{-1}$ min$^{-1}$]');

% plot phase plane       
figure('Color', [1 1 1])
hold on
% shade h >= 0 
fill([patient0.ceMin 20 20 patient0.ceMin], [-10 -10 10 10], [0.9 0.9 0.9], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);            
% shade hDot >= -a1 * h
ceVals = linspace(-20, 20, 500);
bound = -1 * patient0.a1 * (ceVals - patient0.ceMin);
fill([ceVals fliplr(ceVals)], [bound fliplr(ones(size(ceVals)) * 10)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% plot barriers
plot(ceVals, bound, 'k--', 'LineWidth', 1); 
xline(patient0.ceMin, 'k--', 'lower bound for surgery', 'LineWidth', 1);  
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0

% plot phase pictures
xlim([-10 10]); ylim([-10 10]);
plot(xHist1(4,:), ceDotHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(x01(4,1), ceDotHist1(1), 'Color', '#66da54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist2(4,:), ceDotHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(x02(4,1), ceDotHist2(1), 'Color','#cada54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist3(4,:), ceDotHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(x03(4,1), ceDotHist3(1), 'Color','#dcd62b', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist4(4,:), ceDotHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(x04(4,1), ceDotHist4(1), 'Color','#dda010', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist5(4,:), ceDotHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(x05(4,1), ceDotHist5(1), 'Color','#dd5510', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist6(4,:), ceDotHist6, 'LineWidth', 2, 'Color', '#980000');
plot(x06(4,1), ceDotHist6(1), 'Color','#980000', ...
    'Marker', '.', 'MarkerSize', 20);
grid 

legend('1st-order safe set', '2nd-order safe set');
title('Safety phase portrait (exponential controller)')
xlabel('Effect site concentration [$\mu$g mL$^{-1}$]')
ylabel('Effect site concentration time derivative [$\mu$g mL$^{-1}$ min$^{-1}$]')

% 2. simulate graceful controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('RUNNING GRACEFUL CONTROLLER + EKF...\n');
disp(patient0.ceDes);

% print ceG, the failsafe barrier
fprintf('Failsafe concentration barrier ceG [$\mu$g mL$^{-1}$]: %d\n', patient0.ceG);
graceful = true; % use graceful controller 
tStart = 0; tEnd = 60; nTimes = 600;
% scenario 1
disp('Running scenario 1...');
[patient1, tHist, tHistContinuous1, xHist1, xEstHist1,...
    uHist1, uMinHist1, ceDotHist1,...
    BISHist1, BISEstHist1] = patient1.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 2
disp('Running scenario 2...');
[patient2, tHist2, tHistContinuous2, xHist2, xEstHist2,...
    uHist2, uMinHist2, ceDotHist2,...
    BISHist2, BISEstHist2] = patient2.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 3
disp('Running scenario 3...');
[patient3, tHist3, tHistContinuous3, xHist3, xEstHist3,...
    uHist3, uMinHist3, ceDotHist3,...
    BISHist3, BISEstHist3] = patient3.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 4
disp('Running scenario 4...');
[patient4, tHist4, tHistContinuous4, xHist4, xEstHist4,...
    uHist4, uMinHist4, ceDotHist4,...
    BISHist4, BISEstHist4] = patient4.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 5
disp('Running scenario 5...');
[patient5, tHist5, tHistContinuous5, xHist5, xEstHist5,...
    uHist5, uMinHist5, ceDotHist5,...
    BISHist5, BISEstHist5] = patient5.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);

% scenario 6
disp('Running scenario 6...');
[patient6, tHist6, tHistContinuous6, xHist6, xEstHist6,...
    uHist6, uMinHist6, ceDotHist6,...
    BISHist6, BISEstHist6] = patient6.simulate( ...
    graceful, clamp, observer, tStart, tEnd, nTimes);
% 
% 
set(groot, 'defaultAxesTickLabelInterpreter','tex');
% set(groot, 'defaultLegendInterpreter','latex');
% set(0,'defaulttextInterpreter','latex');

% plot effect trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, BISHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, BISHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, BISHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, BISHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, BISHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, BISHist6, 'LineWidth', 2, 'Color', '#980000');
yline(patient0.BISmax, 'k--', 'upper bound for surgery');
yline(patient0.BISdes, 'k--', 'desired');
title('Bispectral index trajectories in time (graceful controller)');
xlabel('Time [min]');
ylabel('Bispectral index [\%]');

% plot site concentration trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, xHist1(4,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xHist2(4,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xHist3(4,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xHist4(4,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xHist5(4,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xHist6(4,:), 'LineWidth', 2, 'Color', '#980000');
yline(patient0.ceG, 'k--', 'padding', 'FontSize', 15);
yline(patient0.ceMin, 'k--', 'lower bound for surgery');
title(['Effect site concentration trajectories in time ' ...
    '(graceful controller)']);
xlabel('Time [min]');
ylabel('Concentration [$\mu$g mL$^{-1}$]');

% plot inputs
figure('Color', [1 1 1])
hold on
plot(tHist, uHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, uHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, uHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, uHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, uHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, uHist6, 'LineWidth', 2, 'Color', '#980000');
yline(uMax, 'k--', 'maximum rate');
title(['Infusion rate trajectories in time ' ...
    '(graceful controller)']);
xlabel('Time [min]');
ylabel('Infusion rate [$\mu$g mL$^{-1}$ min$^{-1}$]');

% plot phase plane        
figure('Color', [1 1 1])
hold on

% shade h >= 0 
fill([patient0.ceMin 100 100 patient0.ceMin], [-50 -50 50 50], [0.7 0.7 0.7], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.3);
% fill([patient0.ceG 100 100 patient0.ceG], [-50 -50 50 50], [0.4 0.4 0.4], ...
%      'EdgeColor', 'none', 'FaceAlpha', 0.6);

% shade new secondary barrier
ceVals = linspace(-100, 100, 500);
lambda1 = -1 * patient0.w * sqrt(patient0.z^2 - 1) + patient0.z*patient0.w; 
lambda2 = patient0.w * sqrt(patient0.z^2 - 1) + patient0.z*patient0.w; 
% set1 = -1 * lambda1 * (ceVals-patient0.ceG)/(patient0.ceG-patient0.ceMin);
% set2 = -1 * lambda2 * (ceVals-patient0.ceG)/(patient0.ceG-patient0.ceMin);
% fill([ceVals fliplr(ceVals)], [set1 fliplr(ones(size(ceVals)) * 10)], ...
%          [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% fill([ceVals fliplr(ceVals)], [set2 fliplr(ones(size(ceVals)) * 10)], ...
%          [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% % plot barriers
% plot(ceVals, set1, 'k--', 'LineWidth', 1); 
% plot(ceVals, set2, 'k--', 'LineWidth', 1); 
xline(patient0.ceMin, 'k--', 'lower bound for surgery', 'LineWidth', 1);
xline(patient0.ceG, 'k--', 'padding', 'LineWidth', 1); 
yline(0, 'k', 'LineWidth', 0.5);  % hDot = 0

% plot 'danger zone' (actuation limits)


% plot phase pictures
xlim([-10 10]); ylim([-10 10]);
plot(xHist1(4,:), ceDotHist1, 'LineWidth', 2, 'Color', '#66da54');
plot(x01(4,1), ceDotHist1(1), 'Color', '#66da54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist2(4,:), ceDotHist2, 'LineWidth', 2, 'Color', '#cada54');
plot(x02(4,1), ceDotHist2(1), 'Color','#cada54', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist3(4,:), ceDotHist3, 'LineWidth', 2, 'Color', '#dcd62b');
plot(x03(4,1), ceDotHist3(1), 'Color','#dcd62b', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist4(4,:), ceDotHist4, 'LineWidth', 2, 'Color', '#dda010');
plot(x04(4,1), ceDotHist4(1), 'Color','#dda010', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist5(4,:), ceDotHist5, 'LineWidth', 2, 'Color', '#dd5510');
plot(x05(4,1), ceDotHist5(1), 'Color','#dd5510', ...
    'Marker', '.', 'MarkerSize', 20);
plot(xHist6(4,:), ceDotHist6, 'LineWidth', 2, 'Color', '#980000');
plot(x06(4,1), ceDotHist6(1), 'Color','#980000', ...
    'Marker', '.', 'MarkerSize', 20);
grid 

legend('Safe set');
title('Safety phase portrait (graceful controller)')
xlabel('Effect site concentration [$\mu$g mL$^{-1}$]')
ylabel('Effect site concentration time derivative [$\mu$g mL$^{-1}$ min$^{-1}$]')


% plot site concentration error trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, xEstHist1(4,:)-xHist1(4,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xEstHist2(4,:)-xHist2(4,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xEstHist3(4,:)-xHist3(4,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xEstHist4(4,:)-xHist4(4,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xEstHist5(4,:)-xHist5(4,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xEstHist6(4,:)-xHist6(4,:), 'LineWidth', 2, 'Color', '#980000');
title('EKF residual trajectories: Effect site concentration');
xlabel('Time [min]');
ylabel('Concentration error [$\mu$g mL$^{-1}$]');

% plot plasma eerror trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, xEstHist1(1,:)-xHist1(1,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xEstHist2(1,:)-xHist2(1,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xEstHist3(1,:)-xHist3(1,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xEstHist4(1,:)-xHist4(1,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xEstHist5(1,:)-xHist5(1,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xEstHist6(1,:)-xHist6(1,:), 'LineWidth', 2, 'Color', '#980000');
title('EKF residual trajectories: Plasma concentration');
xlabel('Time [min]');
ylabel('Concentration error [$\mu$g mL$^{-1}$]');

% plot slow peripheral trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, xEstHist1(3,:)-xHist1(3,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xEstHist2(3,:)-xHist2(3,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xEstHist3(3,:)-xHist3(3,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xEstHist4(3,:)-xHist4(3,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xEstHist5(3,:)-xHist5(3,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xEstHist6(3,:)-xHist6(3,:), 'LineWidth', 2, 'Color', '#980000');
title('EKF residual trajectories: Slow peripherals');
xlabel('Time [min]');
ylabel('Concentration error [$\mu$g mL$^{-1}$]');

% plot fast peripheral trajectories
figure('Color', [1 1 1])
hold on
plot(tHist, xEstHist1(2,:)-xHist1(2,:), 'LineWidth', 2, 'Color', '#66da54');
plot(tHist, xEstHist2(2,:)-xHist2(2,:), 'LineWidth', 2, 'Color', '#cada54');
plot(tHist, xEstHist3(2,:)-xHist3(2,:), 'LineWidth', 2, 'Color', '#dcd62b');
plot(tHist, xEstHist4(2,:)-xHist4(2,:), 'LineWidth', 2, 'Color', '#dda010');
plot(tHist, xEstHist5(2,:)-xHist5(2,:), 'LineWidth', 2, 'Color', '#dd5510');
plot(tHist, xEstHist6(2,:)-xHist6(2,:), 'LineWidth', 2, 'Color', '#980000');
title('EKF residual trajectories: Fast peripherals');
xlabel('Time [min]');
ylabel('Concentration error [$\mu$g mL$^{-1}$]');