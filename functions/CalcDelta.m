function [delta_tot] = CalcDelta(MCMCPar,delta_tot,delta_normX,CR);
% Calculate total normalized Euclidean distance for each crossover value

% Derive sum_p2 for each different CR value 
for zz = 1:MCMCPar.nCR;
    
    % Find which chains are updated with zz/MCMCPar.nCR
    idx = find(CR==zz/MCMCPar.nCR); 
    
    % Add the normalized squared distance tot the current delta_tot;
    delta_tot(1,zz) = delta_tot(1,zz) + sum(delta_normX(idx,1));
    
end;