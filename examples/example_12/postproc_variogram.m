%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                            %
% THIS PROGRAMS HELPS TO POSTPROCESS AND PLOT THE RESULTS OF THE DREAM PACKAGE                               %
%                                                                                                            %
% I would not suggest to use this for the synthetic mathematical problems, but instead only the time series  %
% problems, for instance (simple example) the HYMOD rainfall - runoff example                                %
%                                                                                                            %
% Written by Jasper A. Vrugt                                                                                 %
%                                                                                                            %
% Version 0.5: April 2012: 	Initial setup and evaluation                                                 %
% Version 1.0: May 2012:    Generalization to problems with and without simulation writing, more plotting    %
%                                                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First assemble all sequences in one matrix
ParSet = genparset(Sequences,MCMCPar);

% Find the maximum aposteriori parameter values (last column of ParSet are log-density values!)
idx = find(ParSet(:,end)==max(ParSet(:,end))); idx = idx(1);

% Print those to screen
MAP = ParSet(idx,1:MCMCPar.n)

% Take the last 25% of the posterior samples -- assume that these samples
% are posterior samples (double check that R_stat < 1.2 for all parameters)
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : MCMCPar.n );

% Calculate the mean posterior value of each parameter
MEAN = mean(Pars)

% Calculate the posterior standard deviation of the parameters
STD = std(Pars)

% Calculate the MCMCPar.n-dimensional parameter correlation matrix (R-values)
CORR = corrcoef(Pars)

% Set figure number
fig_number = 1;

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- EVOLUTION OF R_STATISTIC OF GELMAN AND RUBIN ------------------------------
% ------------------------------------------------------------------------------------------------------------

% Now plot the R_statistic for each parameter
figure(fig_number),
% Update figure number
fig_number = fig_number + 1;
% First print the R-statistic of Gelman and Rubin (each parameter a different color)
semilogy(output.R_stat(:,1),output.R_stat(:,2:MCMCPar.n+1)); hold on;
% Add labels
xlabel('Number of iterations','fontsize',14,'fontweight','bold','fontname','Times');
ylabel('R_{stat}','fontsize',14,'fontweight','bold','fontname','Times');
% Add title
title('Convergence of sampled chains','fontsize',14,'fontweight','bold','fontname','Times');
% Now add the theoretical convergence value of 1.2 as horizontal line
plot([0 output.R_stat(end,1)],[1.2 1.2],'k--','linewidth',2);
% Set the the axes
axis([0 output.R_stat(end,1) 1 10]);
% Add a legend
evalstr = strcat('legend(''par.1''');
% Each parameter a different color
for j = 2:MCMCPar.n,
    % Add the other parameters
    evalstr = strcat(evalstr,',''par. ',num2str(j),'''');
end;
% And now conclude with a closing bracket
evalstr = strcat(evalstr,');');
% Now evaluate the legend
eval(evalstr);

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- HISTORGRAMS OF MARGINAL DENSITIES OF PARAMETERS ---------------------------
% ------------------------------------------------------------------------------------------------------------

% Plot the histograms (marginal density) of each parameter;
% What lay out of marginal distributions is desired subplot(r,t)
r = 3; t = 2;
% How many figures do we need to create with this layout?
N_fig = ceil( MCMCPar.n / (r * t) ); counter = 1; j = 1; 
% Open new figure
figure(fig_number);
% Now plot each parameter
while counter <= MCMCPar.n
    % Check whether to open a new figure?
    if j == (r * t) + 1,
        % Update fig_number
        fig_number = fig_number + 1;
        % Open new figure
        figure(fig_number);
        % Reset j to 1
        j = 1;
    end;
    % Now create histogram
    [N,X] = hist(Pars(:,counter));
    % And plot histogram in red
    subplot(r,t,j),bar(X,N/sum(N),'r'); hold on; % --> can be scaled to 1 if using "trapz(X,N)" instead of "sum(N)"!
    if j == 1,
        % Add title
        title('Histograms of marginal distributions of individual parameters','fontsize',14,'fontweight','bold','fontname','Times');
    end;
    % Add x-labels
    evalstr = strcat('Par',{' '},num2str(counter)); xlabel(evalstr,'fontsize',14,'fontweight','bold','fontname','Times');
    % Then add y-label (only if j == 1 or j = r;
    if j == 1 | ( min(abs(j - ([1:r]*t+1))) == 0 ),
        ylabel('Marginal density','fontsize',14,'fontweight','bold','fontname','Times');
    end;
    % Now determine the min and max X values of the plot
    minX = min(X); maxX = max(X); minY = 0; maxY = max(N/sum(N));
    % Now determine appropriate scales
    deltaX = 0.1*(maxX - minX);
    % Calculate x_min and x_max
    x_min = minX - deltaX; x_max = maxX + deltaX;
    % Now determine the min and max Y values of the plot
    y_min = 0; y_max = 1.1*maxY;
    % Lets add the MAP value
    plot(MAP(counter),0.98*y_max,'bx','Markersize',15,'linewidth',3);
    % Adjust the axis
    axis([x_min x_max y_min y_max]);
    % Check if counter = 1,
    if counter == 1, % --> add a title for first figure 
        % Add title
        title('Histograms of marginal distributions of individual parameters','fontsize',14,'fontweight','bold','fontname','Times');
    end;
    % Now update the counter
    counter = counter + 1;
    
    % Update j
    j = j + 1;
end;

% Update fig_number
fig_number = fig_number + 1;

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- CORRELATION PLOTS OF THE POSTERIOR PARAMETER SAMPLES ----------------------
% ------------------------------------------------------------------------------------------------------------

% Open a new plot
figure(fig_number); fig_number = fig_number + 1;
% Plot a matrix (includes unscaled marginals on main diagonal!
plotmatrix(Pars,'+r');
% Add title
title('Marginal distributions and two-dimensional correlation plots of posterior parameter samples','fontsize',14,'fontweight','bold','fontname','Times');

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- CALCULATE THE RMSE OF THE BEST SOLUTION -----------------------------------
% ------------------------------------------------------------------------------------------------------------

% Now compute the RMSE of the best solution, but only if Measurement.MeasData exists!
if exist('Measurement') == 1,

    % Now check whether output simulations have been saved our not?
    if strcmp(MCMCPar.modout,'No');
        % Generate model prediction for best parameter values
        evalstr = ['ModPred = ',ModelName,'(MAP,Extra);']; eval(evalstr);
    else
        % Derive model simulation from fx
        ModPred = fx(idx,1:end);
    end;

    % Compute the RMSE of the maximum aposteriori solution
    RMSE_MAP = sqrt ( sum ( ( ModPred(:) - Measurement.MeasData).^2) / prod(size(ModPred)) )

else

    % Do nothing
    RMSE_MAP = []
end

% If you use option 3, then this RMSE value should be equal to "sqrt(-max(X(:,end-1))/Measurement.N)" !!
% Hence, option 3 uses a standard Gaussian likelihood function, minimizing the SSE (RMSE)

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- CONVERGENCE OF INDIVIDUAL CHAINS TO TARGET DISTRIUBUTION ------------------
% ------------------------------------------------------------------------------------------------------------

% Define colors for different chains
symbol = {'ys','rx','g+','ko','c<'};

% Now loop over each parameter
for j = 1:MCMCPar.n,
    % Open new figures
    figure(fig_number);
    % Update fig_number
    fig_number = fig_number + 1;
    % How many elements does the sequences have
    Nseq = size(Sequences,1) - 1; % --> the first (starting) point is also stored!
    % This gives us the following sample no.
    if MCMCPar.T > 1,
        Sample = [1 MCMCPar.T : MCMCPar.T : MCMCPar.T * Nseq]';
    else
        Sample = [1 : Nseq+1]';
    end;
    % Now plot a number of chains
    for i = 1:min(MCMCPar.seq,5);
        plot(Sample,Sequences(1:end,j,i),char(symbol(i)),'markersize',3,'linewidth',3); if i == 1; hold on; end;
    end
    % Add an axis
    if exist('ParRange'),
        % Use scaling with prior parameter ranges
        axis([0 Nseq + 1 ParRange.minn(j) ParRange.maxn(j)]);
    else
        % Ranges have not been defined -- need to derive them from ParSet
        min_j = min(ParSet(:,j)); max_j = max(ParSet(:,j)); 
        % Now make the ranges a little wider
        if min_j < 0,
            min_j = 1.1*min_j;
        else
            min_j = 0.9*min_j;
        end;
        if max_j > 0,
            max_j = 1.1*max_j;
        else
            max_j = 0.9*max_j;
        end;
        % And scale the figure
        axis([0 Nseq + 1 min_j max_j]);
    end;
    % Lets add the MAP value
    plot(1 * Nseq,MAP(j),'bx','Markersize',15,'linewidth',3);
    % Add a legend
    evalstr = strcat('legend(''chain. 1''');
    % Each parameter a different color
    for jj = 2:min(MCMCPar.seq,5),
        % Add the other parameters
        evalstr = strcat(evalstr,',''chain.',{' '},num2str(jj),'''');
    end;
    % And now conclude with a closing bracket
    evalstr = strcat(evalstr,');');
    % Now evaluate the legend
    eval(char(evalstr));
    % Add a title
    xlabel('Sample number in chain','fontsize',14,'fontweight','bold','fontname','Times');
    % Then add a y-label
    evalstr = strcat('par ',{' '},num2str(j)); ylabel(evalstr,'fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
    % Then add title
    title('Chain convergence plot','fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
end;

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- PLOT THE 95% POSTERIOR SIMULATION UNCERTAINTY -----------------------------
% ------------------------------------------------------------------------------------------------------------

% Calculate the empirical variogram
[H,G,semv,vcloud] = calc_vario(x,y,z);

% How many posterior samples do we like?
N = size(Pars,1);

% Set the separation distance, hd
hd = [0 : 20 : 800];

% Calculate the variogram for each posterior sample
V(1,:) = vmatern(Pars(1,2),Pars(1,3),Pars(1,4),Pars(1,5),hd);

% Now create the remaining part of V
V( 2 : N, 1 : size(hd,2) ) = zeros( N - 1 , size(hd,2) );

% Calculate the posterior variogram ensemble
for i = 2 : N
    % Calculate the variogram for each posterior sample
    V( i , 1 : size(hd,2) ) = vmatern(Pars(i,2),Pars(i,3),Pars(i,4),Pars(i,5),hd);
end

% Calculate the mean posterior variogram
Vmean = mean(V); Vstd = std(V);

% Calculate the median posterior variogram
Vmed = median(V);

% Calculate the 2.5 percentile posterior variogram
Vlo = prctile(V,2.5);

% Calculate the 97.5 percentile posterior variogram
Vhi = prctile(V,97.5);

% Now plot in one figure
figure( fig_number ); plot(hd,Vmean, '-r', hd,Vlo, '--g', hd,Vhi,'--b', semv(:,1),semv(:,2),'ko');
% Add labels
xlabel('Separation distance (lag) [m]','fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
ylabel('Semivariance','fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
title('Posterior variogram prediction uncertainty and empirical variogram','fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
% Add legends
legend('mean variogram','2.5% percentile','97.5% percentile','empirical variogram');