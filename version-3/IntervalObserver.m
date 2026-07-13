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
    
    T = real(T);
    diagD = real(diagD);
    Tinv = inv(T);
    Lz = Tinv * L;
    distBoundZ = abs(Tinv) * (d * ones(4,1));
    noiseBoundZ = abs(Lz) * V;

    % Error in z-world for each state, now with a nonzero initial error:
    % de/dt = lambda * e + w, e(0) = e0 (not necessarily 0)
    % where w = distBoundZ + noiseBoundZ
    % Particular/SS soln: de/dt = 0 --> e_ss = -w/lambda
    % Homogeneous soln: drop w, separate vars --> e_h = C*exp(lambda * t)
    % General soln: e(t) = e_ss + C*exp(lambda * t)
    % Solve for C: e(0) = e0 --> e0 = e_ss + C --> C = e0 - e_ss
    % --> e(t) = e_ss + (e0 - e_ss)*exp(lambda * t)
    %          = e_ss*(1 - exp(lambda * t)) + e0*exp(lambda * t)
    % Worst case (triangle ineq., since e0 and w are magnitude bounds,
    % and exp(lambda*t) > 0 since lambda is real):
    %   width(t) = |e0|*exp(lambda*t) + e_ss*(1 - exp(lambda*t))
    
    e0z = abs(Tinv) * (e0 * ones(4,1)); 
    widthInf = (distBoundZ + noiseBoundZ) ./ (-diagD);
    widthT = e0z .* exp(diagD * tSim) + widthInf .* (1 - exp(diagD * tSim));
    Tpos = max(T, 0); Tneg = min(T, 0);

    % Upper interval width by end
    xHiT = Tpos*widthT + Tneg*(-widthT);
    cost = xHiT(3) + 0.3*sum(xHiT);
end
