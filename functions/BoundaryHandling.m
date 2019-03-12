function [x] = BoundaryHandling(x,ParRange,BoundHandling);
% Function to check whether parameter values remain within prior bounds

% First determine the size of new
[m,n] = size(x);

% Now replicate minn and maxn
minn = repmat(ParRange.minn,m,1); maxn = repmat(ParRange.maxn,m,1);

% Now find which elements of x are smaller than their respective bound
[ii_low] = find(x < minn); 

% Now find which elements of x are larger than their respective bound
[ii_up] = find(x > maxn); 

% Reflection
if strcmp(BoundHandling,'Reflect');

    % reflect in minn
    x(ii_low)= 2 * minn(ii_low) - x(ii_low);     

    % reflect in maxn
    x(ii_up)= 2 * maxn(ii_up) - x(ii_up); 

end;

% Bound
if strcmp(BoundHandling,'Bound');

    % set lower values to minn
    x(ii_low)= minn(ii_low); 
    
    % set upper values to maxn
    x(ii_up)= maxn(ii_up);

end;

% Folding
if strcmp(BoundHandling,'Fold');

    % Fold parameter space lower values
    x(ii_low) = maxn(ii_low) - ( minn(ii_low) - x(ii_low) );
    
    % Fold parameter space upper values
    x(ii_up) = minn(ii_up) + ( x(ii_up) - maxn(ii_up) );

end;

% Now double check in case elements are still out of bound -- this is
% theoretically possible if values are very small or large

% Now double check if all elements are within bounds
[ii_low] = find(x < minn); x(ii_low) = minn(ii_low) + rand(size(ii_low)).* ( maxn(ii_low) - minn(ii_low) );
[ii_up]  = find(x > maxn); x(ii_up)  = minn(ii_up)  + rand(size(ii_up)).*  ( maxn(ii_up)  - minn(ii_up)  );