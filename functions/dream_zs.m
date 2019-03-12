function [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement)

% Calculate MCMCPar.steps
MCMCPar.steps = floor(MCMCPar.steps(MCMCPar)); 

% Calculate MCMCPar.m0
MCMCPar.m0 = MCMCPar.m0(MCMCPar);   

% open an output file with warnings
fid = fopen('warning_file.txt','w');
fprintf(fid,'-------------- DREAM_{(ZS)} warning file --------------\n');

% Check that MCMCPar.m is large enough to create offspring
if MCMCPar.m0 < (2 * MCMCPar.DEpairs * MCMCPar.seq),
    % Warning -- not enough chains to do sampling -- increase number of chains!
    evalstr = char(strcat('DREAM_{(ZS)} WARNING: size of Z not sufficient -> Increase MCMCPar.m0 to at least ',{' '},...
        num2str(2 * MCMCPar.DEpairs * MCMCPar.seq),'\n'));   
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
    % Stop DREAM
    Sequences = []; X = []; output = []; Z = [];
    % And return to main program
    return;
end;
    
% Check how many input variables
if nargin < 5,
    % Define Measurement
    Measurement.MeasData = []; 
end;

if nargin < 4,
    % Specify very large initial parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];
end;

if nargin < 3,
    % Define structure Extra to be empty
    Extra = [];
end;
    
% Calculate the number of calibration data measurements
Measurement.N = size(Measurement.MeasData,1);

% Initialize waitbar
h = waitbar(0,'Running DREAM_{(ZS)} - Please wait...');

% Check whether to restart or not?
if strcmp(MCMCPar.Restart,'No'),

    % Set random number generator
    opts.Seed = 'sum(100*clock)  % evaluated if it is a string';
    % Then generate new seed
    if ischar(opts.Seed)
        randn('state', eval(opts.Seed));     % random number generator state
    else
        randn('state', opts.Seed);
    end

    % Step 0: Initialize variables
    [MCMCPar,pCR,lCR,CR,Iter,iteration,T,fx,m_func,Sequences,Z,Table_JumpRate,iloc,output] = InitMCMC(MCMCPar,Measurement);
    
    % Step 1: Sample MCMCPar.m0 points in the parameter space and store in Z
    
    if strcmp(MCMCPar.prior,'LHS'),
        % Initialize chains with Latin hypercube sampling 
        [Zinit] = LHS(ParRange.minn,ParRange.maxn,MCMCPar.m0 + MCMCPar.seq);
    elseif strcmp(MCMCPar.prior,'COV');
        % Initialize chains with (multi)-normal distribution
        [Zinit] = repmat(MCMCPar.mu,MCMCPar.m0 + MCMCPar.seq,1) + randn(MCMCPar.m0 + MCMCPar.seq,MCMCPar.n) * chol(MCMCPar.cov);
    elseif strcmp(MCMCPar.prior,'PRIOR');
        % Create the initial position of each chain by drawing each parameter individually from the prior
        for qq = 1:MCMCPar.n,
            for zz = 1:MCMCPar.m0 + MCMCPar.seq,
                Zinit(zz,qq) = eval(char(MCMCPar.prior_marginal(qq)));
            end;
        end;
    end;

    % Do boundary handling -- what to do when points fall outside bound
    if strcmp(MCMCPar.BoundHandling,'None') == 0;
        [Zinit] = BoundaryHandling(Zinit,ParRange,MCMCPar.BoundHandling);
    end;

    % Define initial MCMCPar.m0 rows of Z to be initial sample -- posterior density is not needed and thus not evaluated!!
    Z(1:MCMCPar.m0,1:MCMCPar.n) = Zinit(1:MCMCPar.m0,1:MCMCPar.n);

    % Define initial population from last MCMCPar.seq samples of Zinit
    X = Zinit(MCMCPar.m0 + 1:MCMCPar.m0+MCMCPar.seq,1:MCMCPar.n); clear Zinit;

    % Calculate posterior density associated with each value of X
    tic
    [p,log_p,fx(:,1:MCMCPar.seq),MCMCPar.Best] = CompDensity(X,MCMCPar,Measurement,ModelName,Extra);

    % Append X with information about posterior density (or transformation thereof) -- also store model simulations of X
    X = [X p log_p]; Xfx = fx; 

    % Set the first point of each of the MCMCPar.seq sequences equal to the initial X values 
    Sequences(1,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(X',1,MCMCPar.n+2,MCMCPar.seq); 
    
    % Save N_CR in memory
    output.CR(1,1:MCMCPar.nCR+1) = [Iter pCR]; delta_tot = zeros(1,MCMCPar.nCR);

    % Compute the R-statistic of Gelman and Rubin
    [output.R_stat(1,1:MCMCPar.n+1)] = [Iter Gelman(Sequences(1,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];

else % If a restart run is being done: just load the output from the previous ongoing trial

    load DREAM_ZS.mat; 
    %load GWF1P81R04.mat
    MCMCPar.ndraw = 2 * MCMCPar.ndraw;

end;
ETime=toc/3600;                           %Running time
RTime=(MCMCPar.ndraw*ETime/Iter)-ETime;   %Remaining time
disp(['Iter:' num2str(Iter) '  SSE=' num2str(MCMCPar.Best) ...
      '     Time=' num2str(ETime*60) 'min' ...
      '     Elapsed Time=' num2str(ETime) 'hr' '    Remaining Time=' num2str(RTime) 'hr']);
tic;

% Move prior population to posterior population ...
while (Iter < MCMCPar.ndraw),

    % Check that exactly MCMCPar.ndraw are done (uneven numbers this is impossible, but as close as possible)
    if (MCMCPar.steps * MCMCPar.seq) > MCMCPar.ndraw - Iter; 
        % Change MCMCPar.steps in last iteration 
        MCMCPar.steps = ceil((MCMCPar.ndraw - Iter)/MCMCPar.seq);
        % Warning -- not enough chains to do sampling -- increase number of chains!
        evalstr = char(strcat('DREAM_{(ZS)} WARNING: Changed MCMCPar.steps =',{' '},num2str(MCMCPar.steps),{' '},'at',{' '},...
            num2str(Iter),{' '},'function evaluations \n'));
        % Now print warning to screen and to file
        fid = fopen('warning_file.txt','w');
        fprintf(evalstr); fprintf(fid,evalstr);    
    end;
    
    % Initialize totaccept;
    totaccept = 0;

    % Loop a number of times before calculating convergence diagnostic, etc.
    for gen_number = 1:MCMCPar.steps,

        % Update T
        T = T + 1;

        % Define the current locations and associated log-densities
        xold = X(1:MCMCPar.seq,1:MCMCPar.n); log_p_xold = X(1:MCMCPar.seq,MCMCPar.n + 2);

        % Without replacement draw rows from Z for proposal creation
        R = randsample(MCMCPar.m, 2 * MCMCPar.DEpairs * MCMCPar.seq); Zoff = Z(R,1:MCMCPar.n);
        
        % Determine to do parallel direction or snooker update
        if (rand <= MCMCPar.parallelUpdate),
            Update = 'Parallel_Direction_Update';
        else
            Update = 'Snooker_Update';
        end;

        % Generate candidate points (proposal) in each chain using either snooker or parallel direction update
        [xnew,CR(:,gen_number),alfa_s] = offde(xold,Zoff,CR(:,gen_number),MCMCPar,Update,Table_JumpRate,ParRange);

        % Compute the likelihood of each proposal in each chain
        [p_xnew,log_p_xnew,fx_new,MCMCPar.Best] = CompDensity(xnew,MCMCPar,Measurement,ModelName,Extra);

        % Calculate the Metropolis ratio
        [accept] = metrop(MCMCPar,xnew,log_p_xnew,xold,log_p_xold,alfa_s);

        % And update X and the model simulation
        idx_X = find(accept == 1); 
        X(idx_X,1:MCMCPar.n+2) = [xnew(idx_X,1:MCMCPar.n) p_xnew(idx_X) log_p_xnew(idx_X)]; 
        Xfx(:,idx_X) = fx_new(:,idx_X);

        % Check whether to add the current points to the chains or not?
        if (T == MCMCPar.T);
            % Store the current sample in Sequences
            iloc = iloc + 1; Sequences(iloc,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(X',1,MCMCPar.n+2,MCMCPar.seq); 
            % Check whether to store the simulation results of the function evaluations
            if Measurement.N > 0,
                if strcmp(MCMCPar.modout,'Yes'),
                    fx(:,m_func + 1 : m_func + MCMCPar.seq) = Xfx(:,1:MCMCPar.seq);
                    % Update m_func
                    m_func = m_func + MCMCPar.seq;
                end;
            end;
            % And set the T to 0
            T = 0;
        end

        % Compute squared jumping distance for each CR value
        if strcmp(MCMCPar.pCR,'Yes');
            % Calculate the standard deviation of each dimension of X
            r = repmat(std(X(1:MCMCPar.seq,1:MCMCPar.n)),MCMCPar.seq,1);
            % Compute the Euclidean distance between new X and old X
            delta_normX = sum(((xold(1:end,1:MCMCPar.n) - X(1:end,1:MCMCPar.n))./r).^2,2);
            % Use this information to update sum_p2 to update N_CR
            [delta_tot] = CalcDelta(MCMCPar,delta_tot,delta_normX,CR(1:MCMCPar.seq,gen_number));
        end;

        % Check whether to append X to Z
        if (mod(gen_number,MCMCPar.k) == 0),
            % Append X to Z
            Z(MCMCPar.m + 1 : MCMCPar.m + MCMCPar.seq,1:MCMCPar.n+2) = X(1:MCMCPar.seq,1:MCMCPar.n+2);
            % Update MCMCPar.m
            MCMCPar.m = MCMCPar.m + MCMCPar.seq;
        end;

        % How many candidate points have been accepted -- for Acceptance Rate
        totaccept = totaccept + sum(accept);

        % Update Iteration
        Iter = Iter + MCMCPar.seq;

        % Update the waitbar
        waitbar(Iter/MCMCPar.ndraw,h);
        
    end;

    % Reduce MCMCPar.steps to get rounded iteration numbers
    if (iteration == 2), MCMCPar.steps = MCMCPar.steps + 1; end;

    % Store Important Diagnostic information -- Acceptance Rate
    output.AR(iteration,1:2) = [Iter 100 * totaccept/(MCMCPar.steps * MCMCPar.seq)];

    % Store Important Diagnostic information -- Probability of individual crossover values
    output.CR(iteration,1:MCMCPar.nCR+1) = [Iter pCR];

    % Check whether to update individual pCR values
    if (Iter <= 0.1 * MCMCPar.ndraw);
        if strcmp(MCMCPar.pCR,'Yes');
            % Update pCR values
            [pCR] = AdaptpCR(MCMCPar,delta_tot,lCR);
        end;
    end;

    % Generate CR values based on current pCR values
    [CR,lCRnew] = GenCR(MCMCPar,pCR); lCR = lCR + lCRnew;

    % Calculate Gelman and Rubin Convergence Diagnostic
    start_idx = max(1,floor(0.5*iloc)); end_idx = iloc;
    
    % Compute the R-statistic using 50% burn-in from Sequences
    [output.R_stat(iteration,1:MCMCPar.n+1)] = [Iter Gelman(Sequences(start_idx:end_idx,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];

    % Update the iteration
    iteration = iteration + 1;
    Time=toc/3600;                            %Loop time
    ETime=ETime+Time;                         %Elapsed time
    RTime=(MCMCPar.ndraw*ETime/Iter)-ETime;   %Remaining time
    disp(['Iter:' num2str(Iter) '  SSE=' num2str(MCMCPar.Best) ...
        '     Time=' num2str(Time*60) 'min' ...
        '     Elapsed Time=' num2str(ETime) 'hr' '    Remaining Time=' num2str(RTime) 'hr']);
    tic;
    
    % Check whether to save the ouput?
    if strcmp(MCMCPar.save,'Yes');    
    
        % Initialize waitbar
        close(h); h = waitbar(Iter/MCMCPar.ndraw,'Saving iterations to DREAM_{(ZS)}.mat - Please wait...');

        % Store in memory
        save DREAM_ZS.mat
    
        % Rename the waitbar
        close(h); h = waitbar(Iter/MCMCPar.ndraw,'Running DREAM_{(ZS)} - Please wait...');
    end;

end;

% ---------------------------- POST PROCESSING ----------------------------
% Variables have been pre-allocated in memory --> might need remove zeros

% Start with CR
output.CR = output.CR(1:iteration-1,1:size(pCR,2)+1);
% Then R_stat
output.R_stat = output.R_stat(1:iteration-1,1:MCMCPar.n+1);
% Then AR
output.AR = output.AR(1:iteration-1,1:2);
% Then Sequences
Sequences = Sequences(1:iloc,1:MCMCPar.n+2,1:MCMCPar.seq);
% Then the history Z
Z = Z(1:MCMCPar.m,1:MCMCPar.n+2);
% Close the waitbar
close(h);
% Write final line of warning file
fid = fopen('warning_file.txt','w');
fprintf(fid,'----------- End of DREAM_{(ZS)} warning file ----------\n');
% Close the warning file
fclose('all');
% Open the warning file
edit warning_file.txt

% If five output arguments are requested then return fx
if nargout == 5,
    if strcmp(MCMCPar.modout,'Yes'),
       % remove zeros and tranpose fx --> easier to use
       fx = fx(:,1:m_func)';
    else
       fx = [];
    end; 
end;

% -------------------------------------------------------------------------