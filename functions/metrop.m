function [accept] = metrop(MCMCPar,x,log_p_x,xold,log_p_xold,alfa_s);
% Metropolis rule for acceptance or rejection

% And initialize accept with zeros
accept = zeros(MCMCPar.seq,1);

% Calculate the Metropolis ratio
alfa = exp(log_p_x - log_p_xold);

% -------------------------------------------------------------------------
% Using a prior or not? If so, alfa needs to adjusted: multiplied with prior
% -------------------------------------------------------------------------
if strcmp(MCMCPar.prior,'PRIOR');
	
    % Compute prior densities for each parameter in each sequence
    for qq = 1:MCMCPar.n,
        for zz = 1:MCMCPar.seq,
            % Compute prior of proposal 
            prior_x(zz,qq) = eval(char(strrep(MCMCPar.prior_marginal(qq),'rnd(','pdf(x(zz,qq),')));
            % Compute prior of current location
            prior_old(zz,qq) = eval(char(strrep(MCMCPar.prior_marginal(qq),'rnd(','pdf(xold(zz,qq),')));
        end;
    end;

    % Take the product of the various densities
    prior_old = max(prod(prior_old,2),1e-299); prior_prop = max(prod(prior_x,2),1e-299); % (to avoid 0/0 --> NaN)
    % Take the ratio
    alfa_pr = prior_prop./prior_old;
    % Now update alfa value with prior
    alfa = alfa.*alfa_pr;  

end;
% -------------------------------------------------------------------------

% Modify 
alfa = alfa.* alfa_s;

% Generate random numbers
Z = rand(MCMCPar.seq,1);

% Find which alfa's are greater than Z
idx = find(alfa > Z);

% And indicate that these chains have been accepted
accept(idx,1) = 1;