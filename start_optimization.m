% define problem
M = 3; % number of measurements
N = 6; % number of variables

% A{1}: equiangular tight frame
A{1} = [1 1/sqrt(5) 1/sqrt(5) 1/sqrt(5) 1/sqrt(5) 1/sqrt(5);...
        0 2/sqrt(5) 1/10*(5-sqrt(5)) 1/10*(5-sqrt(5)) 1/10*(-5-sqrt(5)) 1/10*(-5-sqrt(5));...
        0 0 sqrt(1/10*(5+sqrt(5))) -sqrt(1/10*(5+sqrt(5))) sqrt(1/10*(5-sqrt(5))) -sqrt(1/10*(5-sqrt(5)))];

% A{2}: subsampled orthogonal matrix
A{2} = [   -0.0386   -0.9101   -0.3310   -0.2553    0.6817    0.7285;...
    0.5994    0.2155    0.9309   -0.4273    0.7302   -0.1092;...
   -0.7995   -0.3539    0.1545   -0.8673   -0.0448   -0.6763];

% A{3}: rotated subsampled orthogonal matrix
A{3} = [ -0.6107   -0.5108    0.0688   -0.0246   -0.5509    0.2394;...
    0.5211   -0.0579   -0.3326    0.4112   -0.3325    0.5787;...
    0.1702    0.2523   -0.0675   -0.8637   -0.2643    0.2950];

% A{4}: realization of random gaussian i.i.d. matrix 
A{4} = [   -0.4087   -0.4408    0.0801    0.7851   -0.1055    0.0685;...
   -0.5769    0.5678    0.1264   -0.2794    0.1539   -0.4766;...
   -0.0015    0.7435   -0.6537    0.1305   -0.0288    0.0443];

p0 = ones(N,1); % x \sim uniform( B_p(1) )

% optimize structured mmse 
% set parameter grid
mvals = 1:4; % number of matrices
dvals = [5,9]; % degree of polynomial nonlinearity, pick odd number
pvals = [0.2:0.2:2];

L = cell(1,3);
[L{:}] = ndgrid(dvals,pvals,mvals);

L = [L{1}(:),L{2}(:),L{3}(:)];
matlabpool open 36
tic
for l = 1:size(L,1)
    [Wcost{l},W{l},aopt{l}] = nonlinear_mmse_altopt(A{L(l,3)},M,N,L(l,1),L(l,2)*p0, zeros(L(l,1)+1,1), 10*pinv(A{L(l,3)}) );
end
etime = toc
matlabpool close


% optimize linear mmse
% set parameter grid
mvals = 1:4; % number of matrices
dvals = 1; % number of iterations
pvals = [0.2:0.2:2];

L0 = cell(1,3);
[L0{:}] = ndgrid(dvals,pvals,mvals);

L0 = [L0{1}(:),L0{2}(:),L0{3}(:)];
for l = 1:size(L0,1)
    [Wcost0{l},W0{l},aopt0{l}] = nonlinear_mmse_altopt(A{L0(l,3)},M,N,L0(l,1),L0(l,2)*p0,[0;1], 10*pinv(A{L0(l,3)}) );
end

save data_aws