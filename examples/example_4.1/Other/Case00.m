clear all 
clc

%Original case ISOR,3HC,1VHC,1BC,3PR,1RCH = 10 variables 
xi=[1;NaN;5e-8;NaN;0;10;10;-20;5E-5];

%[1] Transmisivity
for n=1:4
    if n==2
        HC=[ones(9,5)*7e-2,[ones(4,4)*1e-2;ones(5,4)*1e-2]];
        dlmwrite(eval(['''Tran' num2str(n,'%1i') '.dat''']), ones(9,9).*HC, 'delimiter', '\t', 'precision', 7);
    elseif n==4
        tic
        x=(0.5:1:9)*500;                            % Spatial coordinates of the cell centers
        y=(0.5:1:9)*500;                            % Spatial coordinates of the cell centers
        [X,Y]=meshgrid(x,y);
        a=X(:);
        b=Y(:);
        HC=zeros(size(a,1),1);
        for m=1:size(a,1)
            HC(m,1)=-20*pi*cos(pi*a(m))*sin(pi*b(m))...
                -20*pi*sin(pi*a(m))*cos(pi*b(m))...
                +40*pi^2*(1+a(m)+b(m))*cos(pi*a(m))*cos(pi*b(m));
        end
        HC=reshape(HC,9,9)*1e-8; 
        dlmwrite(eval(['''Tran' num2str(n,'%1i') '.dat''']), ones(9,9).*HC, 'delimiter', '\t', 'precision', 7);
        toc
    else
        dlmwrite(eval(['''Tran' num2str(n,'%1i') '.dat''']), ones(9,9)*xi(n,1), 'delimiter', '\t', 'precision', 7);
    end
end


%[2] Boundary condition
for n=5:5
    dlmwrite(eval(['''BC' num2str(n-4,'%1i') '.dat''']), [ones(9,1)*xi(n,1) zeros(9,8)], 'delimiter', '\t', 'precision', 7);
end
close('all')

%[3] Wells
wel=fopen('twri.wel','w');
fprintf(wel,'         3         0    MXWELL,IWELBD\n');
fprintf(wel,'         3              ITMP (NWELLS)\n');
fprintf(wel,'         1         3        8       %6.2f\n',xi(6,1));
fprintf(wel,'         1         4        8       %6.2f\n',xi(7,1));
fprintf(wel,'         1         5        2       %6.2f\n',xi(8,1));
fclose(wel);

%[4] Recharge
rch=fopen('twri.rch','w');
fprintf(rch,'         1      0      NRCHOP,IRCHBD\n');
fprintf(rch,'         1             INRECH\n');
fprintf(rch,'         0%7.3E       RECH-1\n',xi(9,1));
fclose(rch);

%[5] MODFLOW
!mf2005.exe twri.nam

   
%[6]Plotting results    
[X,Y]=meshgrid((0.5:1:9)*500,(0.5:1:9)*500);  % Cell centers
Hfid=fopen('head.hed','r');
TrueHead=fscanf(Hfid,'%f\n');
fclose(Hfid);
TrueHead=reshape(TrueHead,9,18)';

SP(:,:,1)=TrueHead(1:9,:);
SP(:,:,2)=textread('Tran2.dat');
SP(:,:,3)=TrueHead(10:end,:);
SP(:,:,4)=textread('Tran4.dat');
figure(1)
for n=1:4
subplot(2,2,n)
surf(X,Y,SP(:,:,n))
xlim([min(X(:)) max(X(:))])
ylim([min(Y(:)) max(Y(:))])
end
min(min(SP(:,:,4)))
save TrueHead.mat TrueHead 


%[7] Save Files for Tecplot
dx=500;
dy=500;
nx=9;
ny=9;
[X,Y]=meshgrid((0:nx)*dx,(0:ny)*dy);
for n=1:2
    HC=SP(:,:,n);   %[1] HC [2] Head
    HC=[HC HC(:,end);HC(end,:) HC(end,end)]; %#ok<AGROW>
    D1=[X(:) Y(:) HC(:)];
    if n==1
        fid1=fopen(eval(['''HC' num2str(n,'%1i') '.dat''']),'wt');
    elseif n==2
        fid1=fopen(eval(['''HEAD' num2str(n-1,'%1i') '.dat''']),'wt');
    end
    fprintf(fid1,'Title="K1"\n');
    fprintf(fid1,'Variables = X Y V\n');
    fprintf(fid1, 'Zone T="K1" i = %1.0f , j = %1.0f,  k = %1.0f\n',10,10,1);
    fprintf(fid1, '%10.0f %10.0f %7.5f\n', D1');
    fclose(fid1);
end


    