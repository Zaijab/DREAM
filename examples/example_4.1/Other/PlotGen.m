function PlotGen(x,gen,fitness)
clf
rangemin=-3;
rangemax=3;
OptRange=[rangemin;rangemax];
ParaRange=[1e-4;1e-1];    xi=NaN(length(x),1);
for m=1:length(x)
    xi(m,1) = interp1(OptRange,ParaRange,x(m,1));
end
HC=zeros(9,9,1);
HC(:,:,1)=reshape(xi(1:81),9,9);
HC(:,:,2)=reshape(xi(82:end),9,9);
[X,Y]=meshgrid((0.5:1:9)*500,(0.5:1:9)*500);  % Cell centers
figure(1)
for n=1:2
    subplot(1,2,n)
    surf(X,Y,HC(:,:,n))
    xlim([min(X(:)) max(X(:))])
    ylim([min(Y(:)) max(Y(:))])
    zlim([0 0.08])
    title(eval(['''Generation:' num2str(gen) '   Fitness:' num2str(fitness) '''']))
end
drawnow
end