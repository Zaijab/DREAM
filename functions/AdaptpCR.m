function [pCR] = AdaptpCR(MCMCPar,delta_tot,lCR); 
% Updates the probabilities of the various crossover values

% Adapt pCR using information from averaged normalized jumping distance
pCR = MCMCPar.seq * (delta_tot./lCR) / sum(delta_tot);

% Normalize pCR
pCR = pCR./sum(pCR);