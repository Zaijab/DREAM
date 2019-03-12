% Postprocessing of DREAM results and visualization of results

% Plot R-statistic
figure
hold on
for i = 1:size(Sequences,2)-2
    plot(output.R_stat(:,1),output.R_stat(:,i+1))
end
ind = find(isfinite(output.R_stat(:,end)));
plot(output.R_stat([ind(1) end],1),[MCMCPar.Rstatcon MCMCPar.Rstatcon],'r')


% Generate ParSet discarding the first ndisc draws (burn-in as indicated by R-statistic plot)
nburnin = 1000;
start = floor(nburnin/MCMCPar.seq);
stop = size(Sequences,1);

ParSet = zeros((stop-start+1)*size(Sequences,3),size(Sequences,2)-1);
for i = 1:size(Sequences,3)
	ParSet((i-1)*(stop-start+1)+1:i*(stop-start+1),:) = Sequences(start:stop,1:size(Sequences,2)-1,i);
end

% Plot histograms of posterior parameter sample
figure
for i = 1:length(ParRange.minn)
	subplot(length(ParRange.minn),1,i)
	hist(ParSet(:,i))
	xlim([ParRange.minn(i) ParRange.maxn(i)])
end
