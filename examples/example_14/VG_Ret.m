function [RMSE,p_pred] = VG_Ret(pars,x,p);
% Now predict function

alpha = pars(1); n = pars(2);

% predict
p_pred = (1 + (alpha * x).^n).^(1/n - 1);

% Calculate RMSE
RMSE = sqrt ( sum ( ( p_pred - p).^2) / ( prod(size(x)) ) );