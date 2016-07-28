% numerical results for optimized performance
vander = @(x,D) ( x(:)*ones(1,D+1) ).^(ones(length(x),1)*(0:D));

% load optimized smmse estimators
load data_aws.mat
niter = 1e5;

for l = 1:size(L0,1)
    X{l} = gen_vec(N,L0(l,2),niter,0);
    Xhat{l} = zeros(N,niter);
    for i = 1:niter
        Xhat{l}(:,i) = vander(W0{l}*A{L0(l,3)}*X{l}(:,i),L0(l,1))*aopt0{l};
    end
    disp(l)
end

for l = 1:size(L0,1)
    lmmse_num(l) = sum(sum((X{l}-Xhat{l}).^2))/niter;
end

save('data_lmmse.mat','lmmse_num')