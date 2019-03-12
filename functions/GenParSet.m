function ParSet = GenParSet(Sequences,MCMCPar);
% Generates a 2D matrix ParSet from 3D array Sequences

% Determine how many elements in Sequences
[NrX,NrY,nrZ] = size(Sequences); 

% Initalize ParSet
ParSet = [];

% If save in memory -> No -- ParSet is empty
if (NrX == 0),
    % Do nothing
else
    % ParSet derived from all sequences
    for qq = 1:MCMCPar.seq,
        ParSet = [ParSet; Sequences(:,:,qq) (1:size(Sequences(:,:,qq),1))'];
    end;
    ParSet = sortrows(ParSet,[MCMCPar.n+3]); ParSet = ParSet(:,1:MCMCPar.n+2);
end;