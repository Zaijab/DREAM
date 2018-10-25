function [ER1,ER2,xn] = excess(x_loss,cmax,bexp,Pval,PETval);
% this function calculates excess precipitation and evaporation

xn_prev = x_loss;
ct_prev = cmax*(1-power((1-((bexp+1)*(xn_prev)/cmax)),(1/(bexp+1))));
% Calculate Effective rainfall 1
ER1 = max((Pval-cmax+ct_prev),0.0);
Pval = Pval-ER1;
dummy = min(((ct_prev+Pval)/cmax),1);
xn = (cmax/(bexp+1))*(1-power((1-dummy),(bexp+1)));
% Calculate Effective rainfall 2
ER2 = max(Pval-(xn-xn_prev),0);

% Alternative approach
evap = (1-(((cmax/(bexp+1))-xn)/(cmax/(bexp+1))))*PETval; % actual ET is linearly related to the soil moisture state
xn = max(xn-evap, 0); % update state

%evap = min(xn,PETval);
%xn = xn-evap;