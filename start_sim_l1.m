% baseline using solutions to l1
% define problem
M = 3; % number of measurements
N = 6; % number of variables

A{1} = [1 1/sqrt(5) 1/sqrt(5) 1/sqrt(5) 1/sqrt(5) 1/sqrt(5);...
        0 2/sqrt(5) 1/10*(5-sqrt(5)) 1/10*(5-sqrt(5)) 1/10*(-5-sqrt(5)) 1/10*(-5-sqrt(5));...
        0 0 sqrt(1/10*(5+sqrt(5))) -sqrt(1/10*(5+sqrt(5))) sqrt(1/10*(5-sqrt(5))) -sqrt(1/10*(5-sqrt(5)))];
    
A{2} = [   -0.0386   -0.9101   -0.3310   -0.2553    0.6817    0.7285;...
    0.5994    0.2155    0.9309   -0.4273    0.7302   -0.1092;...
   -0.7995   -0.3539    0.1545   -0.8673   -0.0448   -0.6763];

A{3} = [ -0.6107   -0.5108    0.0688   -0.0246   -0.5509    0.2394;...
    0.5211   -0.0579   -0.3326    0.4112   -0.3325    0.5787;...
    0.1702    0.2523   -0.0675   -0.8637   -0.2643    0.2950];

A{4} = [   -0.4087   -0.4408    0.0801    0.7851   -0.1055    0.0685;...
   -0.5769    0.5678    0.1264   -0.2794    0.1539   -0.4766;...
   -0.0015    0.7435   -0.6537    0.1305   -0.0288    0.0443];

p0 = ones(N,1); % x \sim uniform( B_p(1) )

% optimize structured mmse 
% set parameter grid
mvals = [1,3,4]; % number of matrices
dvals = 0; 
pvals = [0.2:0.2:1];

niter = 2500;

L2 = cell(1,3);
[L2{:}] = ndgrid(dvals,pvals,mvals);
L2 = [L2{1}(:),L2{2}(:),L2{3}(:)];

warning('off','MATLAB:bitcmp:NDiscontinueSupport');
for l = 1:size(L2,1)
    X{l} = gen_vec(N,L2(l,2),niter,0);
    for i = 1:niter
        if L2(l,2)<=1 
            cvx_begin quiet
               cvx_solver sdpt3
               cvx_precision high
               variable x(N)
               minimize( norm( x, 1 ) )
               subject to
                    A{L2(l,3)} * x == A{L2(l,3)}*X{l}(:,i);
            cvx_end
            Xhat{l}(:,i) = x;
        else
            Xhat{l}(:,i) = pinv(A{L2(l,3)})*A{L2(l,3)}*X{l}(:,i);
        end
    end
    disp(l)
end

for l = 1:size(L2,1)
    mmse_l1(l) = sum(sum((X{l}-Xhat{l}).^2))/niter;
end

save('data_l1.mat','mmse_l1','L2')