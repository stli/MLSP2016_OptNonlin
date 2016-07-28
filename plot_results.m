% run scripts to generate results
plot_true = 1;

% run the following scripts in order to generate the results
% for simulation of start_optimization.m an AWS C4.8XLARGE instance was used
% start_optimization();
% start_sim_smmse();
% start_sim_lmmse();
% start_sim_l1();

% compute signal power E[x.'x] for normalization
M = 3; N = 6;
Px = [];
for p = 0.2:0.2:2
    EX = gen_matrices_EX([],M,N,[],p*ones(N,1),[]);
    Px(end+1) = trace(EX);
end

% load analytic errors for smmse
load('data_aws.mat','Wcost','Wcost0','pvals','L','L0')
for i = 1:length(Wcost)
    smmse_ana(i) = Wcost{i}(end);
end
for i = 1:length(Wcost0)
    lmmse_ana(i) = Wcost0{i}(end);
end

% load numeric errors for smmse
load('data_smmse')

% load numeric errors for lmmse
load('data_lmmse')

% load numeric errors for l1/l2
load('data_l1','mmse_l1','L2')

% % set colors for plots
if plot_true
	try
		addpath('C:\Users\limmer\Documents\matlab\toolbox\cbrewer')
		CT=cbrewer('div','RdYlBu', 4);
		CT(3,:) = CT(2,:);
	catch
		CT = colormap;
	end
    MK = {'o','*','s','v'};


    for m=[1,3,4] % plot results for matrices A_1,...,A_4
        % get indices of vectors to plot
        idx_D9 = find( ismember(L(:,[1,3]),[9,m],'rows' ) );
        idx_lin= find( ismember(L0(:,[3]),[m],'rows' ) );
        idx_l12= find( ismember(L2(:,[3]),[m],'rows' ) );

    %% plot smmse
        plot(pvals(2:end),smmse_num(idx_D9(2:end))./Px(2:end),'-','LineWidth',4,'Linestyle','--','Color',[CT(m,:)]);   
        hold on
        grid on
        axis tight
        h1(m) = plot(pvals(2:end),smmse_ana(idx_D9(2:end))./Px(2:end),horzcat('-',MK{m}),'LineWidth', 2,'Color',[CT(m,:)],...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[CT(m,:)],...        
            'MarkerSize',7);
    %% plot l12
        plot(pvals(2:5),mmse_l1(idx_l12(2:5))./Px(2:5),horzcat(':',MK{m}),'LineWidth', 2,'Color',[CT(m,:)],...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[CT(m,:)],...
            'MarkerSize',7);
    end

      h2 =  plot(pvals(2:end),lmmse_ana(idx_lin(2:end))./Px(2:end),'--','LineWidth', 2,'Color','k',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',3);
    set(0,'DefaultAxesColorOrder',CT)
    set(gca,'TickLabelInterpreter', 'latex')
    set(gca,'DefaultTextInterpreter', 'latex')

    legend([h1(1) h1(4) h1(3) h2],{'equiangular tight frame','normalized iid','subsampled orthogonal','linear estimation bound'},...
        'Location','SouthEast', 'Interpreter', 'latex')
    axis tight
    grid on
    xlabel('$p$', 'Interpreter', 'latex')
    ylabel('normalized MSE')
end