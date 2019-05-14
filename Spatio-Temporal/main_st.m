% @author: Hessam Babaee
clc
clear all
close all


addpath '../Utilities/'
MS = 'MarkerSize';
set(0,'defaulttextinterpreter','latex')

R     = [216, 82,  24 ]/255;
B     = [0 , 113, 188]/255;  

global ModelInfo
load('../Data/Years-2105-2016')
%% Setup
N_L = size(LF.X,1);
N_H = size(HF.X,1);
D   = size(LF.X,2);

jitter = eps;
ModelInfo.jitter=jitter;

%% Generate Data
ModelInfo.X_H = HF.X;
ModelInfo.y_H = HF.Y;

ModelInfo.X_L = LF.X;
ModelInfo.y_L = LF.Y;

%% Training
hyp = [log([1 1 1 1 1 1 1 1]) 1 -1 -1];
options = optimoptions('fminunc','GradObj','on','Display','iter',...
    'Algorithm','trust-region','Diagnostics','on','DerivativeCheck','on',...
    'FinDiffType','central');
[ModelInfo.hyp,~,~,~,~,~] = fminunc(@likelihood,hyp,options);


%% Predictions
load(['../Data/SST-Sattelite-',num2str(Year)])
load('../Data/Station_info')

% LSeg = 10000;
% NSeg = floor(size(Pred(1).X,1)/LSeg);
% n = 0;
% [imax, jmax] =size(Pred(1).X2d);
% xmin = min(min(Pred(1).X2d));
% xmax = max(max(Pred(1).X2d));
% ymin = min(min(Pred(1).Y2d));
% ymax = max(max(Pred(1).Y2d));
% Np  = imax*jmax;
%
%
% mean_f_star = [];
% var_f_star  = [];
% for i=1:NSeg
%     i
%     if (i~=NSeg)
%         [M, V] = predictor_f_H(Pred(1).X((i-1)*LSeg+1:i*LSeg,:));
%     else
%         [M, V] = predictor_f_H(Pred(1).X((i-1)*LSeg+1:end,:));
%     end
%     mean_f_star = [mean_f_star; M];
%     var_f_star = [var_f_star; V];
% end






N = 2*365;
for i=1:Pred(2).NPart
    X = Pred(2).X;
    Y = Pred(2).Y;
    I = Pred(2).Part;
    
    x_s = [X(I(i,1),1)*ones(N,1) X(I(i,1),2)*ones(N,1) [1:N]'+t_ref];
    time = x_s(:,3)-t_ref;
    [M, V] = predictor_f_H(x_s);
    
    
    x_t_H= ConvertTimetoDays(X(I(i,1):I(i,2),3),t_ref);
    y_t_H = Y(I(i,1):I(i,2),:);
    
    
    X = Pred(3).X;
    Y = Pred(3).Y;
    I = Pred(3).Part;
    
    x_t_L= ConvertTimetoDays(X(I(i,1):I(i,2),3),t_ref);
    y_t_L = Y(I(i,1):I(i,2),:);
    
    
    
    %% Plot Results
    clear h;
    h = figure
    clear leg;
    hold
    if (ismember(n_trial(i),n_train))
        leg1 = 'MWRA (Training)';
    else
        leg1 = 'MWRA (Test)';
    end
    c1 = B;
    
    leg2 = 'Satellite (Test)';
    
    h(1) = plot(x_t_H, y_t_H,'o','color','k',MS,10,'MarkerFaceColor', 'k');
    h(2) = plot(x_t_L, y_t_L,'Square','color',R,MS,10);
    h(3) = plot(time,M,'LineWidth',3,'color',B);
    [l,h(4)] = boundedline(time, M, 2.0*sqrt(V), ':', 'alpha','cmap', c1);
    outlinebounds(l,h(4));
    
    
    leg{1} = leg1;
    leg{2} = leg2;
    leg{3} = sprintf('Multifidelity GP');
    leg{4} = sprintf('Two standard deviation');
    
    hl = legend(h,leg,'Location','southeast');
    legend boxoff
    set(hl,'Interpreter','latex')
    xlabel('Days')
    ylabel('SST ($^{\circ}$ C)')
    ylim([-5 25])
    grid
    
    title(['Station:',Station.S(n_trial(i))'])
    axis square
    xlim([1 2*365]);
    set(gca,'FontSize',15);
    set(gcf, 'Color', 'w');
    
    
end

%% --------- Validation: NERACOOS Buoy A01 -----------
figure
lat = 42.52;
lon= -70.56;
[A1,txt,raw] = xlsread('../Data/A1-2015.xlsx');
N_A1 = size(A1,1);
X_A1 = [repmat([lon lat], N_A1,1) A1(:,1)];
[M_A1, V_A1] = predictor_f_H(X_A1);

[B1,txt,raw] = xlsread('../Data/A1-2016.xlsx');
N_B1 = size(B1,1);
X_B1 = [repmat([lon lat], N_B1,1) B1(:,1)];
[M_B1, V_B1] = predictor_f_H(X_B1);

A1 = [A1; B1];
M_A1=[M_A1; M_B1];
V_A1=[V_A1; V_B1];

plot(A1(:,1)-t_ref, A1(:,2),'o','color',R,MS,6);
hold on
plot(A1(:,1)-t_ref,M_A1,'LineWidth',3,'color',B);
boundedline(A1(:,1)-t_ref, M_A1, 2.0*sqrt(V_A1), ':', 'alpha','cmap', c1);

xlabel('Days')
ylabel('SST ($^{\circ}$ C)')
ylim([-5 25])
grid


axis square
xlim([1 2*365]);
set(gca,'FontSize',15);
set(gcf, 'Color', 'w');

hl = legend('NERACOOS Buoy A01','Multifidelity GP','Two standard deviation','Location','southeast')
legend boxoff
set(hl,'Interpreter','latex')


t=[1:10:300]'+t_ref;
X = [repmat([lon lat], length(t),1) t];
[M, V] = predictor_f_H(X);
M_A1=[];
for i=1:length(t)
    [m,I]=min(abs(A1(:,1)-t(i,1)));
    M_A1(i) = A1(I,2);
end


rmpath ../Utilities