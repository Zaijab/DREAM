function [H,G,semv,vcloud]=calc_vario(xd,yd,z,nlag,hmax,tolag)
% calculate distance matrix & empirical variogram
%   Input:
%   xd,yd   : x & y coordinate
%   z       : data
%   nlag (optional) : no of lags
%   hmax (optional) : max distance for variogram calculation
%   tolag (optional)   : percent of lag tolerance, from 0 to 100, 
%                       bigger value creates smoother variogram
%   Output:
%   H       : distance matrix, for all points
%   G       : semivariance matrix for all points
%   semv    : empirical semivariogram
%               [distance, semivariance, no. points, std.dev of semivar]
%   vcloud  : variogram cloud
%               [distance, semivariance]
%   Budiman (01.2004)
if(nargin<4),
    nlag=30;
    hmax=0;
    tolag=0;
elseif(nargin<5),
    hmax=0;
    tolag=0;
elseif(nargin<6),
    tolag=0;
end

X=[xd,yd];
[n, nx] = size(X); % n=number of rows in distance matrix

% Based on Peter Acklam distance matrix calculation without loop
% calculate distance matrix
x = permute(X, [ 1 3 2 ]);
y = permute(X, [ 3 1 2 ]);
H = sqrt(sum((x(:,ones(1,n),:) - y(ones(1,n),:,:)).^2, 3));

% calculate semivariance matrix
z1 = permute(z, [ 1 3 2 ]);
z2 = permute(z, [ 3 1 2 ]);
G = 0.5.*(sum((z1(:,ones(1,n),:) - z2(ones(1,n),:,:)).^2, 3));

% extract dist & gamma & make into 1 column, to get variogram cloud
ij=find(tril(G)>0);
hd=H(ij);
gd=G(ij);
vcloud=[hd,gd];

% calculate empirical variogram, method of moments
if(hmax==0),hmax=0.8*max(hd);end;
step=hmax/(nlag+1);     %   step length
xtol=step*tolag/100;    %   lag tolerance
h=[0:step:hmax];
nlag=length(h);
semv=zeros(nlag,4);
for i=1:nlag-1,
    ll=h(i)-xtol;
    ul=h(i)+step+xtol;
    ij=find(hd>=ll & hd<=ul);
    nij=length(ij);
    if(nij>0),
        semv(i,1)=mean(hd(ij));     % mean distance
        semv(i,2)=mean(gd(ij));     % mean semivariance
        semv(i,3)=std(gd(ij));      % std. dev of the semivariance
        semv(i,4)=nij;              % no. of pairs
    end
end

ij=find(semv(:,4)>1);
semv=semv(ij,:);
