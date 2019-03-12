function [DEversion] = DEStrategy(MCMCPar);
% Determine which sequences to evolve with what DE strategy

% Determine probability of selecting a given number of pairs;
p_pair = (1/MCMCPar.DEpairs) * ones(1,MCMCPar.DEpairs); p_pair = cumsum(p_pair); p_pair = [0 p_pair];
% Generate a random number between 0 and 1
Z = rand(MCMCPar.seq,1);
% Select number of pairs
for qq = 1:MCMCPar.seq,
    z = find(Z(qq,1)>p_pair); DEversion(qq,1) = z(end);
end;