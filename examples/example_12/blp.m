function [z0,s0]=blp(x0,y0,x,y,z,c0,c1,range,nu,M,beta,m,CI)
% best linear predictor
% z0 = predicted
% s0 = variance of predicted
X=[x,y];
% distance between unknown & data
X0=[x0,y0];
h=distmat(X, X0);
c=cova_matern(c0,c1,range,nu,h);

% predicted
p1=c'*CI*(z-M*beta);
p2=m'*beta;
z0(:,1)=p1+p2;
z0(:,2)=p1;
z0(:,3)=p2;

% variance of predicted
g=m-M'*CI*c;
e1=(c0+c1)-c'*CI*c;
e2=g'*inv(M'*CI*M)*g;
s0=e1+e2;
s0(:,2)=e1;
s0(:,3)=e2;