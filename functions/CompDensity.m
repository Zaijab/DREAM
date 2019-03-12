function [p,log_p,fx,Best] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra);
% This function computes the density of each x value

% Check whether to store the output of each model evaluation (function call)
if Measurement.N > 0,
    if strcmp(MCMCPar.modout,'Yes'),
        % Create initial ModPred
        fx = NaN(Measurement.N,size(x,1));
    end;
end;

% Loop over the individual parameter combinations of x
for ii = 1:MCMCPar.seq,

    % Call model to generate simulated data (Simulation results for the
    % Likihood function calculations
    evalstr = ['fx(:,ii) = ',ModelName,'(x(ii,:),Extra);']; 
    eval(evalstr);

    % If we have measured data --> calculate the residual (used by most likelihood functions)
    if Measurement.N > 0,
        % Calculate the error residual
        Err = (Measurement.MeasData(:) - fx(1:Measurement.N,ii));
    end;

    if MCMCPar.lik == 1, % Model returns posterior density
        p(ii,1) = fx(1,ii); log_p(ii,1) = log(p(ii,1));
    end;

    if MCMCPar.lik == 2, % Log-density function

        % Derive the log density
        if size(Measurement.Sigma,1) == 1,  % --> homoscedastic error
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - Measurement.N * log( Measurement.Sigma ) - 1/2 * Measurement.Sigma^(-2) * sum ( Err.^2 );
        else                                % --> heteroscedastic error
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - sum ( log( Measurement.Sigma ) ) - 1/2 * sum ( ( Err./Measurement.Sigma ).^2);
        end;

        % And retain in memory
        p(ii,1) = log_p(ii,1);

    end;

    if MCMCPar.lik == 3, % Model returns vector of predictions

        % Derive the sum of squared error
        SSR = sum(abs(Err).^2);
        MCMCPar.Best(MCMCPar.Best>SSR)=SSR;
        Best= MCMCPar.Best;

        % And retain in memory
        p(ii,1) = -SSR; log_p(ii,1) = - Measurement.N/2 * log(SSR);

    end;

    if MCMCPar.lik == 4, % Model returns log posterior density

        p(ii,1) = fx(1,ii); log_p(ii,1) = p(ii,1);

    end;

    if MCMCPar.lik == 7, % Log likelihood with AR-1 model of residuals

        % First order autoregressive (AR-1) correction of residuals
        rho = x(ii,MCMCPar.n); Err_2 = Err(2:Measurement.N,1) - rho * Err(1:Measurement.N-1,1);

        % Now compute the log-likelihood
        if size(Measurement.Sigma,1) == 1,  % --> homoscedastic error
            log_p(ii,1) = -(Measurement.N/2) * log(2*pi) - (Measurement.N/2) * log(Measurement.Sigma^2 / (1-rho^2)) - (1/2) * (1-rho^2) * (Err(1)/Measurement.Sigma)^2 ...
                - (1/2) * sum((Err_2./Measurement.Sigma).^2);
        else                                % --> heteroscedastic error
            log_p(ii,1) = -(Measurement.N/2) * log(2*pi) - (Measurement.N/2) * log(mean(Measurement.Sigma.^2) / (1-rho^2)) - (1/2) * (1-rho^2) * (Err(1)/Measurement.Sigma(1))^2 ...
                - (1/2) * sum((Err_2./Measurement.Sigma(2:Measurement.N)).^2);
        end;

        % And retain in memory
        p(ii,1) = log_p(ii,1);
    end;

    if MCMCPar.lik == 8, % Generalized log likelihood (GL)
        % Extract statistical model parameters
        par = Extra.fpar;               % fixed parameters
        par(Extra.idx_vpar) = x(ii,:);  % variable parameters
        par = par';                     % make it a column vector
        statpar = par(end-10:end);

        % Compute the log-likelihood
        log_p(ii,1) = GL('est',statpar,fx(1:Measurement.N,ii),Measurement.MeasData);

        % And retain in memory
        p(ii,1) = log_p(ii,1);
    end;

    if MCMCPar.lik == 9, % Whittle likelihood function

        % Calculate the log-likelihood using spectral density
        [log_L] = Whittle_logL(Measurement,fx(1:Measurement.N,ii));

        % Now store in memory
        p(ii,1) = log_L; log_p(ii,1) = p(ii,1);

    end

    if MCMCPar.lik == 10, % Approximate Bayesian Computation

        % Now calculate rho
        rho = MCMCPar.rho( fx(1:Measurement.N,ii) , Measurement.MeasData(:) ) + normrnd(0,MCMCPar.delta);

        % How many elements does rho consist of?
        N_rho = prod(size(rho));

        % Easier to work with log-density in practice --> when distance to
        % 0 is large (with small delta)
        p(ii,1) = - ( N_rho / 2) * log(2 * pi) - N_rho * log( MCMCPar.delta ) - 1/2 * MCMCPar.delta^(-2) * sum ( rho.^2 );

        % And log density
        log_p(ii,1) = p(ii,1);

    end;

end;