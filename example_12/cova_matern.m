function S=cova_matern(c0,c1,range,nu,H)
% covariance of Matern

Hr=H./range;
if(nu==0.5),    % run exponential if nu = 0.5
    S=c1.*exp(-Hr);
else
    r1=(2^(nu-1))*gamma(nu);
    bes=besselk(nu,Hr);
    F=(1./r1).*(Hr.^nu).*bes;
    S=c1.*F;
end
i=find(H==0);
S(i)=c0+c1;