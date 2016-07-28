% numerical results for optimized performance
vander = @(x,D) ( x(:)*ones(1,D+1) ).^(ones(length(x),1)*(0:D));

% load optimized smmse estimators
load data_aws.mat
niter = 1e5;

for l = 1:size(L,1)
    X{l} = gen_vec(N,L(l,2),niter,0);
    Xhat{l} = zeros(N,niter);
    for i = 1:niter
        Xhat{l}(:,i) = vander(W{l}*A{L(l,3)}*X{l}(:,i),L(l,1))*aopt{l};
    end
    disp(l)
end

for l = 1:size(L,1)
    smmse_num(l) = sum(sum((X{l}-Xhat{l}).^2))/niter;
end

save('data_smmse.mat','smmse_num')