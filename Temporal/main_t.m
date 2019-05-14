% @author: Hessam Babaee
clc
clear all
close all

addpath ../Utilities
MS = 'MarkerSize';
set(0,'defaulttextinterpreter','latex')
R     = [216, 82,  24 ]/255; 
B     = [0 , 113, 188]/255;  

global ModelInfo
load('../Data/Time-StationF02')
% load('../Data/Time-StationF23')
% load('../Data/Time-StationN18')
%% Setup
N_L = size(LF.X,1);
N_H = size(HF.X,1);
D   = size(LF.X,2);

jitter = eps;
ModelInfo.jitter=jitter;

%% Data
ModelInfo.X_H = HF.X;
ModelInfo.y_H = HF.Y;

ModelInfo.X_L = LF.X;
ModelInfo.y_L = LF.Y;

%% Training
hyp = [log([1 1 1 1]) 1 -1 -1];
options = optimoptions('fminunc','GradObj','on','Display','iter',...
    'Algorithm','trust-region','Diagnostics','on','DerivativeCheck','on',...
    'FinDiffType','central');
[ModelInfo.hyp,~,~,~,~,~] = fminunc(@likelihood,hyp,options);


%% Predictions
load('../Data/Station_info')

x_s = Pred(1).X;
time = x_s(:,end);
[mean_f_star, var_f_star] = predictor_f_H(x_s);

x_t = Pred(2).X;
y_t = Pred(2).Y;

%% Plot Results
clear h;
clear leg;
hold

leg1 = 'MWRA (Training)';
leg2 = 'Satellite (Training)';
h(1) = plot(HF.X-t_ref, HF.Y,'*','color',B,MS,10);
h(2) = plot(LF.X-t_ref, LF.Y,'Square','color',R,MS,10);
h(3) = plot(x_t-t_ref, y_t,'o',MS,10,'color','k','MarkerFaceColor', 'k');
h(4) = plot(x_s-t_ref,mean_f_star,'LineWidth',3,'color',B);
[l,h(5)] = boundedline(x_s-t_ref, mean_f_star, 2.0*sqrt(var_f_star), ':', 'alpha','cmap', B);
outlinebounds(l,h(5));


leg{1} = leg1;
leg{2} = leg2;
leg{3} = 'MWRA (Test)';
leg{4} = sprintf('Multifidelity GP');
leg{5} = sprintf('Two standard deviation');

hl = legend(h,leg,'Location','southeast');
legend boxoff
set(hl,'Interpreter','latex')
xlabel('Days')
ylabel('SST ($^{\circ}$ C)')
grid
axis square
xlim([1 2*365]);
ylim([-5 25]);
set(gca,'FontSize',15);
set(gcf, 'Color', 'w');

rmpath ../Utilities