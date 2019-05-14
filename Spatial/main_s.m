% @author: Hessam Babaee
clc
clear all
close all

addpath ../Utilities
MS = 'MarkerSize';
set(0,'defaulttextinterpreter','latex')
R     = [216, 82,  24 ]/255; 
B     = [0 , 113, 188]/255;  
stateName = 'Massachusetts';
states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
ma = states(strcmp(states.Name, stateName));


global ModelInfo
% load('../Data/May18-2016')
% load('../Data/March20-2015')
load('../Data/March20-2015-withF29')
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
hyp = [log([1 1 1 1 1 1]) 1 -1 -1];
options = optimoptions('fminunc','GradObj','on','Display','iter',...
    'Algorithm','trust-region','Diagnostics','on','DerivativeCheck','on',...
    'FinDiffType','central');
[ModelInfo.hyp,~,~,~,~,~] = fminunc(@likelihood,hyp,options);

%% Predictions
load(['../Data/SST-Sattelite-',num2str(Year)])
load('../Data/Station_info')


[mean_f_star,var_f_star]= predictor_f_H(Pred(1).X);



[imax, jmax] =size(Pred(1).X2d);
Np = imax*jmax;

xmin = min(min(Pred(1).X2d));
xmax = max(max(Pred(1).X2d));
ymin = min(min(Pred(1).Y2d));
ymax = max(max(Pred(1).Y2d));


X = reshape(Pred(1).X(:,1),imax,jmax);
Y = reshape(Pred(1).X(:,2),imax,jmax);

I_L = inpolygon(X(:),Y(:),ma.Longitude,ma.Latitude);
I_C = find(X(:)<-70.6 & Y(:)<41.8);

T_m = mean_f_star;
T_V = sqrt(var_f_star);

T_m(find(I_L==1))=nan;
T_m(I_C)=nan;
T_m = reshape(T_m,imax, jmax);


T_V(find(I_L==1))=nan;
T_V(I_C) = nan;
T_V = reshape(T_V,imax, jmax);

i=1;
for n=n_train
    Stn.S(i)   = Station.S(n);
    Stn.lat(i) = Station.lat(n);
    Stn.lon(i) = Station.lon(n);
    i = i +1;
end 

%% ------------------------------------------------------------------------
figure
contourf(X,Y,T_m,21); view([0 0 1]); axis equal; axis tight; hold on
colorbar('southoutside')
plot(ma.Longitude,ma.Latitude,'LineWidth',2); hold off
axis([xmin xmax ymin ymax])
xlabel('Longitude'); ylabel('Latitude')
set(gca,'FontSize',15);
grid
title('Multi-Fidelity Model: Mean')
%% ------------------------------------------------------------------------
figure
T = Data.SST(10:32,1:25,Day);
surf(Data.X(10:32,1:25),Data.Y(10:32,1:25),T); view([0 0 1]); axis equal; axis tight; hold on
colorbar('southoutside')
plot(ma.Longitude,ma.Latitude,'LineWidth',2);
axis([xmin xmax ymin ymax])
xlabel('Longitude'); ylabel('Latitude')
set(gca,'FontSize',15);
grid
title('Low-Fidelity Data: Satellite')
%% ------------------------------------------------------------------------
figure
contourf(X,Y,T_V,20);  axis equal; axis tight; hold on
colorbar('southoutside')
plot(Stn.lon,Stn.lat,'o',...
    'MarkerFaceColor',R,...
    'MarkerSize',7); 
text(Stn.lon+.01,Stn.lat,Stn.S,'FontSize',15);
alpha(.5)
plot(ma.Longitude,ma.Latitude,'LineWidth',2,'color','k'); hold off
axis([xmin xmax ymin ymax])
xlabel('Longitude'); ylabel('Latitude')
set(gca,'FontSize',15);
grid
title('Multi-Fidelity Model: Uncertainty Map')

rmpath ../Utilities