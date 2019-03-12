function [ii] = PickProp(MCMCPar,log_p_y);
% Pick the proposal based on the weights

% Determine weights of each individual proposal using log-density value
for zz = 1:MCMCPar.seq,
    % Find the appropriate indices
    i_start = (zz-1) * MCMCPar.k_pc + 1; i_end = zz * MCMCPar.k_pc;
    % Scale log density with maximum to avoid problems with zero 
    m = max(log_p_y(i_start:i_end,1)); 
    % Determine the weight --> density
    w(1:MCMCPar.k_pc,zz) = exp(log_p_y(i_start:i_end,1) - m);
    % Renormalize the weights
    w(:,zz) = w(:,zz)./sum(w(:,zz));
end;
% Take the cumulative sum
w = [zeros(1,MCMCPar.seq) ; cumsum(w)];

% Generate a uniform random number between 0 and 1 in each chain
u = rand(1,MCMCPar.seq);

% Match u and w to select appropriate proposal point in each chain
for zz = 1:MCMCPar.seq,    
    % Get i
    idx = find(u(1,zz) > w(:,zz)); idx = idx(end) + (zz-1) * MCMCPar.k_pc;
    % Save idx
    ii(zz,1) = idx;
end; 