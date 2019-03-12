dd=x(2)-x(1);
nn=length(Extra.datasim);

% Check convergence statistics
figure(1)
numb=nnz(output.R_stat(:,1));
plot(output.R_stat(1:numb,1),output.R_stat(1:numb,2:MCMCPar.n));
hold on
plot(output.R_stat(1:numb,1),1.2*ones(numb,1),'r','LineWidth',3);
hold off
title('G-R convergence stat.')
ylabel('Gelman-Rubin R');
xlabel('Function evaluation');
ylim([0.8 3]);

figure(2)
numb=nnz(output.AR(:,1));
plot(output.AR(1:numb,1),output.AR(1:numb,2),'k');
title('Acceptance rate')
ylabel('Acceptance rate (%)');
xlabel('Function evaluation');

% plot true model
figure(3)
model=zeros(30,30);
for j=1:nz
    for i=1:nz
        model(j,i)=slowtrue((j-1)*nz+i);
    end
end
model=1./model;
for k=1:Extra.dimver
    zpatch(1,1)=z(k);
    zpatch(2,1)=z(k)+dd;
    zpatch(3,1)=z(k)+dd;
    zpatch(4,1)=z(k);
    for i=1:Extra.dimhor
        % Option to plot proposed water content models
        xpatch(1,1)=x(i);
        xpatch(2,1)=x(i);
        xpatch(3,1)=x(i)+dd;
        xpatch(4,1)=x(i)+dd;
        patch(xpatch(:,1),-zpatch(:,1),1000*model(k,i),'LineStyle','none')
    end
end

xlabel('Distance (m)')
ylabel('Depth (m)')
ylim([-3 0])
axis('equal')
title('True model (m/micro.s)')
caxis([50 170])
colorbar

for ii=4:2:8
    titles=num2str(ii);
    titles=strcat('Realization order: ', titles);
    figure(ii)
    for kk=1:3
        for jj=1:3
            % plot example model and true model
            model=zeros(30,30);
            for j=1:ii
                for i=1:ii
                    model(j,i)=Sequences(round(numb*(0.6+0.2*(jj-1))),(j-1)*Extra.parx+i,kk);
                end
            end
            model=idct2(model);
            model=10.^model; % Slowness field
            model=1./model;
            %model=(0.3*model-Extra.por*sqrt(Extra.pa)-(1-Extra.por)*sqrt(Extra.ps))/(sqrt(Extra.pw)-sqrt(Extra.pa)); % Water content
            subplot(3,3,(kk-1)*3+jj)
            for k=1:Extra.dimver
                zpatch(1,1)=z(k);
                zpatch(2,1)=z(k)+dd;
                zpatch(3,1)=z(k)+dd;
                zpatch(4,1)=z(k);
                for i=1:Extra.dimhor
                    % Option to plot proposed water content models
                    xpatch(1,1)=x(i);
                    xpatch(2,1)=x(i);
                    xpatch(3,1)=x(i)+dd;
                    xpatch(4,1)=x(i)+dd;
                    patch(xpatch(:,1),-zpatch(:,1),1000*model(k,i),'LineStyle','none')
                end
            end

            xlabel('Distance (m)')
            ylabel('Depth (m)')
            ylim([-3 0])
            axis('equal')
           
            title(titles)
            caxis([50 170])
            colorbar
        end
    end
end
