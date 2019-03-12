function V=vmatern(c0,c1,range,nu,H)
% variogram Matern
% Input:
%   H       : distance

Hr=H./range;
r1=2.^(nu-1)*gamma(nu);
bes=besselk(nu,Hr);
C=(1./r1).*(Hr.^nu).*bes;

V=c0+c1.*(1-C);
ij=find(H==0);
if(length(ij)>0),V(ij)=c0;end;