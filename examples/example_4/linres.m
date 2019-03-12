function [x_slow,outflow] = linres(x_slow,inflow,Rs);
% Linear reservoir
x_slow = (1-Rs)*x_slow + (1-Rs)*inflow;
outflow = (Rs/(1-Rs))*x_slow;
