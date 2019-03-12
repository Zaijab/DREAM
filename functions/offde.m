function [xnew,CR,alfa_s] = offde(xold,Zoff,CR,MCMCPar,Update,Table_JumpRate,ParRange)
% Calculate candidate points using discrete proposal distribution

% Determine how many pairs to use for each jump in each chain
[DEversion] = DEStrategy(MCMCPar);

% Generate uniform random numbers for each chain to determine which dimension to update
D = rand(MCMCPar.seq,MCMCPar.n);

% Generate noise to ensure ergodicity for each individual chain
noise_x = MCMCPar.eps * (2 * rand(MCMCPar.seq,MCMCPar.n) - 1);

% Initialize the delta update to zero
delta_x = zeros(MCMCPar.seq,MCMCPar.n);

if strcmp(Update,'Parallel_Direction_Update'), % PARALLEL DIRECTION UPDATE

    % Define which points of Zoff to use to generate jumps
    rr(1,1) = 1; rr(1,2) = rr(1,1) + DEversion(1) - 1; rr(1,3) = rr(1,2) + 1; rr(1,4) = rr(1,3) + DEversion(1) - 1;
    % Do this for each chain
    for qq = 2:MCMCPar.seq,
        % Define rr to be used for population evolution
        rr(qq,1) = rr(qq-1,4) + 1; rr(qq,2) = rr(qq,1) + DEversion(qq,1) - 1; rr(qq,3) = rr(qq,2) + 1; rr(qq,4) = ...
            rr(qq,3) + DEversion(qq,1) - 1;
    end;

    % Each chain evolves using information from other chains to create offspring
    for qq = 1:MCMCPar.seq,

        % ------------ WHICH DIMENSIONS TO UPDATE? USE CROSSOVER ----------
        [i] = find(D(qq,1:MCMCPar.n) > (1-CR(qq,1)));

        % Update at least one dimension
        if isempty(i), i = randperm(MCMCPar.n); i = i(1); end;
        % -----------------------------------------------------------------

        % Select the appropriate JumpRate and create a jump
        if ( rand < (1 - MCMCPar.pJumpRate_one) ),

            % Now determine gamma, the jump factor
            switch MCMCPar.ABC

                case {'No'}

                    % Select the JumpRate (dependent of NrDim and number of pairs)
                    NrDim = size(i,2); JumpRate = Table_JumpRate(NrDim,DEversion(qq,1));

                case {'Yes'}

                    % Turner (2012) paper -- CU[0.5,1]
                    JumpRate = 0.5 + rand/2;

            end;

            % Produce the difference of the pairs used for population evolution
            delta = sum(Zoff(rr(qq,1):rr(qq,2),1:MCMCPar.n) - Zoff(rr(qq,3):rr(qq,4),1:MCMCPar.n),1);
            % Then fill update the dimension
            delta_x(qq,i) = (1 + noise_x(qq,i)) * JumpRate.*delta(1,i);
        else
            % Set the JumpRate to 1 and overwrite CR and DEversion
            JumpRate = 1; CR(qq,1) = -1;
            % Compute delta from one pair
            delta = Zoff(rr(qq,1),1:MCMCPar.n) - Zoff(rr(qq,4),1:MCMCPar.n);
            % Now jumprate to facilitate jumping from one mode to the other in all dimensions
            delta_x(qq,1:MCMCPar.n) = JumpRate * delta;
        end;
    end;
end;

if strcmp(Update,'Snooker_Update'), % SNOOKER UPDATE
    % Determine the number of rows of Zoff
    NZoff = size(Zoff,1);
    % Define rr
    rr = [1:1:NZoff]; rr = reshape(rr,2,size(rr,2)/2)';
    % Define JumpRate -- uniform rand number between 1.2 and 2.2
    Gamma = 1.2 + rand;
    % Loop over the individual chains
    for qq = 1:MCMCPar.seq,
        % Define which points of Zoff z_r1, z_r2
        zR1 = Zoff(rr(qq,1),1:MCMCPar.n); zR2 = Zoff(rr(qq,2),1:MCMCPar.n);
        % Now select z from Zoff; z cannot be zR1 and zR2
        r = [1:NZoff]; r(rr(qq,1)) = 0; r(rr(qq,2)) = 0; r = r(r>0); t = randperm(NZoff-2);
        % Define z
        z(qq,1:MCMCPar.n) = Zoff(r(t(1)),1:MCMCPar.n);
        % Define projection vector x(qq) - z
        F = xold(qq,1:MCMCPar.n) - z(qq,1:MCMCPar.n); D = max(F*F',1e-300);
        % Orthogonally project of zR1 and zR2 onto F
        zP = F * (sum((zR1-zR2).*F) / D);
        % And define the jump
        delta_x(qq,1:MCMCPar.n) = Gamma * zP;
        % Update CR because we only consider full dimensional updates
        CR(qq,1) = 1;
    end;
end;

% Now propose new x
xnew = xold + delta_x;

% Define alfa_s
if strcmp(Update,'Snooker_Update'),
    % Determine Euclidean distance
    alfa_s = [(sum((xnew - z).^2,2)./sum((xold - z).^2,2)).^((MCMCPar.n-1)/2)];
else
    alfa_s = ones(MCMCPar.seq,1);
end;

% Do boundary handling -- what to do when points fall outside bound
if strcmp(MCMCPar.BoundHandling,'None') == 0;
    [xnew] = BoundaryHandling(xnew,ParRange,MCMCPar.BoundHandling);
end;