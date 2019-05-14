% @author: Maziar Raissi

function [mean_star, var_star] = predictor_f_H(X_star)

global ModelInfo

% X_L = ModelInfo.X_L;
% X_H = ModelInfo.X_H;
% y_L = ModelInfo.y_L;
% y_H = ModelInfo.y_H;

X_L = Normalize(ModelInfo.X_L,min(ModelInfo.X_L),max(ModelInfo.X_L),mean(ModelInfo.X_L));
X_H = Normalize(ModelInfo.X_H,min(ModelInfo.X_L),max(ModelInfo.X_L),mean(ModelInfo.X_L));
y_L = Normalize(ModelInfo.y_L,min(ModelInfo.y_L),max(ModelInfo.y_L),mean(ModelInfo.y_L));
y_H = Normalize(ModelInfo.y_H,min(ModelInfo.y_L),max(ModelInfo.y_L),mean(ModelInfo.y_L));


x_star = Normalize(X_star,min(ModelInfo.X_L),max(ModelInfo.X_L),mean(ModelInfo.X_L));
% x_star = X_star;

hyp = ModelInfo.hyp;
rho = hyp(end-2);

D = size(X_H,2);

y = [y_L; y_H];

L=ModelInfo.L;

psi1 = rho*k(x_star, X_L, hyp(1:D+1),0);
psi2 = rho^2*k(x_star, X_H, hyp(1:D+1),0) + k(x_star, X_H, hyp(D+2:2*D+2),0);
psi = [psi1 psi2];

% calculate prediction
mean_star = psi*(L'\(L\y));

var_star = rho^2*k(x_star, x_star, hyp(1:D+1),0) ...
  + k(x_star, x_star, hyp(D+2:2*D+2),0) ...
  - psi*(L'\(L\psi'));

var_star = abs(diag(var_star));

mean_star = DeNormalize(mean_star,min(ModelInfo.y_L),  max(ModelInfo.y_L),mean(ModelInfo.y_L));
var_star =    var_star *(max(ModelInfo.y_L) - min(ModelInfo.y_L)).^2;
