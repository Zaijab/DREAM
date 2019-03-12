%regression_station
ilon     = 2;
ilat     = 3;
idepth   = 4;
itemp    = 5;
iCN      = 18;
iabsdepth= 19;
istation = 20;


modelFun = @(p,x) p(1)*exp(p(2)*x);
beta = [7 0.001];

clear Result;
Result = zeros(max(CN(:,istation)),12);


for i=1:max(CN(:,istation))
    clear ind;
    ind = find(CN(:,istation) == i & CN(:,iCN) > 2 & CN(:,iCN) < 20 & CN(:,itemp) > -5 & CN(:,itemp) < 33 & CN(:,ilat) >= -90);
    if length(ind) > 3
        X = CN(ind,idepth);
        y = CN(ind,iCN);
        clear paramEsts;
        paramEsts = NonLinearModel.fit(X, y, modelFun, beta);
        Result(i,1) = double(paramEsts.Coefficients(1,1));%intercept
        Result(i,2) = double(paramEsts.Coefficients(2,1));%slope
        Result(i,3) = paramEsts.RMSE;%RMSE
        Result(i,4) = double(paramEsts.Coefficients(1,4));%p value intercept
        Result(i,5) = double(paramEsts.Coefficients(2,4));%p value slope
        X_lin = [ones(size(X)) X];
        [b,bint,r,rint,stats] = regress(y,X_lin);
        Result(i,6) = b(1);
        Result(i,7) = b(2);
        Result(i,8) = stats(3);
        Result(i,9) = CN(ind(1),ilat);%lat
        Result(i,10)= CN(ind(1),ilon);%lon
        
        y_T = CN(ind,itemp);
        bT = regress(y_T,X_lin);
        Result(i,11)= bT(1);%predicted surface T
        Result(i,12) = CN(ind(1),iabsdepth);%abs depth
    else
        Result(i,1:12) = [-9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999  -9999];
    end
end