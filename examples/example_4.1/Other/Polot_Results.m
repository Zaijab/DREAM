clear all
clc
load OD_C1
[obj row]=min(outdat(:,1));
for m=1:10:1500
    %[2] Unscaling CMAES parameters to  MODFLOW parameters
    x=outdat(m,2:end)';
    rangemin=-2;
    rangemax=2;
    OptRange=[rangemin;rangemax];
    ParaRange=[1e-3;1e-1];    xi=NaN(length(x),1);
    for m=1:length(x)
        xi(m,1) = interp1(OptRange,ParaRange,x(m,1));
    end
    HC=reshape(xi,9,9);
    
    %[3] Transmisivity for Layer 1
    fidT2=fopen('Tran2.dat','w');
    for n=1:9
        fprintf(fidT2,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n',HC(n,:));
    end
    fclose(fidT2);
    !mf2005.exe twri.nam > screenout
    
    [X,Y]=meshgrid((0.5:1:9)*500,(0.5:1:9)*500);  % Cell centers
    
    Hfid=fopen('head.hed','r');
    SimuHead=fscanf(Hfid,'%f\n');
    fclose(Hfid);
    SimuHead=reshape(SimuHead,9,18)';
    
    figure(2)
    surf(X,Y,HC)
    xlim([min(X(:)) max(X(:))])
    ylim([min(Y(:)) max(Y(:))])
        
end
load TrueHead
f=sqrt(sum(sum((TrueHead-SimuHead).^2)))













