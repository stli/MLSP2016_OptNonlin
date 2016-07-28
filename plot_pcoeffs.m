% load analytic errors for smmse
load('data_aws.mat','aopt','A','W','pvals','L','N')

lmax = [];
idx_D9 = find( L(:,1) == 9); % get indices for D=9
% find sup norm(W*A*x,'inf') numerically
for i = 1:length(idx_D9)
    Aopt(:,i) = aopt{idx_D9(i)}(:);
    X = gen_vec(N,L(idx_D9(i),2),5e3,0);
    Y = (W{idx_D9(i)}*A{L(idx_D9(i),3)})*X;
    lmax(end+1) = max(abs(Y(:)));
end

% set colors for plots
try
	addpath('C:\Users\limmer\Documents\matlab\toolbox\cbrewer')
	CT=cbrewer('div','RdYlBu', 4);
	colormap(CT)
catch
	CT = colormap;
end
%plot nonlinearity
t = sym('t')
figure
%L(idx_D9(2:5),:)
for i = 2:5
    h(i) = ezplot(poly2sym(flipud(Aopt(:,i)).'),[-lmax(i),lmax(i)])
    title([])
    hold on
end
set(0,'DefaultAxesColorOrder',CT)
set(gca,'TickLabelInterpreter', 'latex')
set(gca,'DefaultTextInterpreter', 'latex')
grid on
legend([h(2:5)],{'$p=0.4$','$p=0.6$','$p=0.8$','$p=1$'},...
    'Location','SouthEast', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\mathcal{T}(t)$')
set(h(2),'Linewidth',3.5)
set(h(3),'Linewidth',3)
set(h(4),'Linewidth',2.5)
set(h(5),'Linewidth',2)
axis([-4 4 -1 1])

