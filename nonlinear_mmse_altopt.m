function [Wcost,W,aopt] = nonlinear_mmse_altopt(A,M,N,D,p,aopt,W)
% alternating optimization for nonlinear mmse optimization as outlined in
% $addreference$
% inputs:
% A = sening matrix
% M = number of measurements
% N = number of variables x_n, n \in \{1,...,N\}
% D = max degree of polynomial nonlinearity
% p = characteristic vector of B_p
% aopt = initial vector of polynomial coefficients
% W0 = initial linear recovery mapping
%
% outputs:
% Wcost = 

% import manopt for the optimization of W
addpath(genpath('C:\Users\limmer\Documents\matlab\toolbox\manopt'))
addpath(genpath('/home/centos/Desktop/matlab/manopt'))

% Create the problem structure.
manifold = euclideanfactory(N,M);
problem.M = manifold;

% Define the problem cost function and its gradient.

% precompute matrices of monomial exponents for faster optimization
for d = 0:2*D
    if d == 0
        alpha{d+1} = zeros(1,N);
    else
        [alpha{d+1}] = multinom_exp(N,d,'descend');
    end
end

% define problem cost
problem.cost = @mmse_cost;
function [f, store] = mmse_cost(W, store)
    
        [EX,EX1,EX2] = gen_matrices_EX(W*A,M,N,D,p,alpha);
        f = trace(EX) - 2*ones(1,N)*EX1*aopt + aopt.'*EX2*aopt;
end

% define problem gradient w.r.t. W
problem.egrad = @mmse_grad;
function [g, store] = mmse_grad(W, store )

    [EXW1,EXW2] = gen_matrices_EXW(W*A,M,N,D,p,alpha,aopt);
    g = (-2*EXW1 + EXW2)*A.';
end
    

% Numerically check gradient and Hessian consistency.
%figure;
%checkgradient(problem);

% Solve.
opts.maxiter = 1e3;
opts.tolgradnorm = 1e-9;
opts.tolcost = 1e-6;
opts.minstepsize = 1e-6;
opts.verbosity = 0;
maxaltiter = 25;

% initialize W0
Wcost(1) = inf;

for i = 1:maxaltiter
    % optimize for a
    [EX,EX1,EX2] = gen_matrices_EX(W*A,M,N,D,p,alpha);
    aopt = pinv(EX2)*(EX1.')*ones(N,1);
    Wcost(2*i) = trace(EX) - 2*ones(1,N)*EX1*aopt + aopt.'*EX2*aopt;
    %disp(['icost: ',num2str(Wcost(end))])
    
    % optimize for W
    [W, Wcost(2*i+1), info] = steepestdescent(problem,W,opts);          %#ok<ASGLU>
    %disp(['icost: ',num2str(Wcost(end))])

    % check convergence
    costdiff1 = Wcost(2*i-1) - Wcost(2*i); costdiff2 = Wcost(2*i) - Wcost(2*i+1);
    if (costdiff1 <= 1e-12) && (costdiff2 <=1e-12)
        break
    end
end
disp(['cost: ',num2str(Wcost(end))])

end