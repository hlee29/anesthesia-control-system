%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntervalObserver.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One-time computation of observer gains for PropofolSystem.m 
% using pole placement in diagonalized basis, 
% where poles are chosen using PSO to minimize interval width
% at the end of the simulation time. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Parameters and observability check

   % Deposition volumes [mL]
    V1 = 6.28 * 1000;               % plasma 
    V2 = 25.5 * 1000;               % fast periph.
    V3 = 273 * 1000;                % slow periph.

    % Clearances [mL⁻¹ min⁻¹]
    CL = 1.79 * 1000;               % plasma
    Q2 = 1.75 * 1000;               % fast periph.
    Q3 = 1.11 * 1000;               % slow periph.

    % Transport coefficients [min⁻¹]
    k10 = CL/V1;                    % plasma -> elimination
    k12 = Q2/V1;                    % fast periph. -> plasma
    k13 = Q3/V1;                    % slow periph. -> plasma
    k21 = Q2/V2;                    % plasma -> fast periph.
    k31 = Q3/V3;                    % plasma -> slow periph. 
    ke0 = 0.146;                    % plasma -> effect site

    A = [-k10-k12-k13   k12     k13     0
         k21          -k21     0       0
         k31           0      -k31     0
         ke0           0       0       -ke0];
    B = [1/V1; 0; 0; 0;]; 
    C = [0 0 0 1]; 
    L = [0; 0; 0; 0.6227]; 
    

    O = [C; C*A;C*A^2;C*A^3];
    disp("RANK OF OBSERVABILITY MATRIX: ")
    rank(O)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. PSO to search for poles that would minimize interval width by the 
% end of th esimulation

d = 0.05;                       % disturbance magnitude bound [μg mL⁻¹]
V = 0.5;                        % sensor noise magnitude bound [μg mL⁻¹]
e0 = 0.3;                       % initial estimation error bound [μg mL⁻¹]
                                % (matches initError in runPropofolSystem.m)
tSim = 70;                      % simulation horizon [min]
                                % tEnd in runPropofolSystem.m

lb = [-5, -5, -5, -5];
ub = [-1e-4, -1e-4, -1e-4, -1e-4];
psOptions = optimoptions('particleswarm');
psOptions.UseParallel = false;  
                               
psOptions.SwarmSize = 100;
psOptions.MaxIterations = 150;
psOptions.MaxStallIterations = 30;
costFun = @(p) poleBoundCost(p, A, C, d, V, e0, tSim);
bestPoles = sort(particleswarm(costFun, 4, lb, ub, psOptions));

L = place(A',C',bestPoles)';
L

A_LC = A - L*C;
eig(A_LC)

[T, D] = eig(A_LC)
Tinv = inv(T)
Tinv*A_LC*T

% Cost function for PSO

function cost = poleBoundCost(poles, A, C, d, V, e0, tSim)
    poles = sort(poles(:)');

    % Enforce distinct eigenvalues
    if min(diff(poles)) < 1e-4
        cost = 1e6;
        return
    end

    % Place the poles
    try
        L = place(A', C', poles)';
    catch
        cost = 1e6;
        return
    end

    A_LC = A - L*C;
    [T, D] = eig(A_LC);
    diagD = diag(D);

    % Check for real and negative eigenvalues (Hurwitz)
    if any(abs(imag(diagD)) > 1e-8) || any(real(diagD) >= 0)
        cost = 1e6;
        return
    end

    % Let interval width W := e_hi + e_lo. Small w denotes its components. 
    % Then, in the z-coordinates, dWz/dt = DWz + 2(|Tinv*L|*V + d_hi - d_lo)
    % Denote by E the last term, 2(Tinv*|L|*V + d_hi - d_lo). It is bound by 2*(distBoundZ + noiseBoundZ).
    % Solving the ODE, we get a steady-state solution wz_ss = -E/lambda,
    % and the homogeneous solution wz_h = C*exp(lambda*t). 
    % So Wz(t) = Wz_ss(t) + Wz_h(t).
    % We have the initial condition Wz(0) = 2e0z, so C = 2e0z - Wz_ss(t) = 2e0z + E/lambda
    % For the final solution Wz(t) = (2e0z + E/lambda)*exp(lambda*t) - E/lambda = 
    % = 2e0z * exp(lambda*t) + Wz_ss(1 - exp(lambda*t)).  
    % Back to x-coordinates, Wx(t) = |T|Wz(t). 

    T = real(T);
    diagD = real(diagD);
    Tinv = inv(T);
    Lz = Tinv * L;
    distBoundZ = abs(Tinv) * (d * ones(4,1));
    noiseBoundZ = abs(Lz) * V;
    
    e0z = abs(Tinv) * (e0 * ones(4,1)); 
    Wz_ss = 2 * (distBoundZ + noiseBoundZ) ./ (-diagD);
    Wz_T = 2 * e0z .* exp(diagD * tSim) + Wz_ss .* (1 - exp(diagD * tSim));

    % Transform back
    Wx_T = abs(T) * Wz_T
    cost = Wx_T(3) + 0.3*sum(Wx_T); % Emphasize slow third state
end
