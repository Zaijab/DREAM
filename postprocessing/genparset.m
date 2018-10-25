function ParSet = GenParSet(Sequences,MCMCPar);
% Generates ParSet

[NrX,NrY,nrZ] = size(Sequences); ParSet = [];

% If save in memory -> No -- ParSet is empty
if (NrX == 1),
    % Do nothing
else
    % If save in memory -> Yes -- ParSet derived from all sequences
    for qq = 1:MCMCPar.seq,
        ParSet = [ParSet; Sequences(:,:,qq) (1:size(Sequences(:,:,qq),1))'];
    end;
    ParSet = sortrows(ParSet,[MCMCPar.n+3]); ParSet = ParSet(:,1:MCMCPar.n+2);
end;