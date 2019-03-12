function f=Case02(x)

fclose('all');
%[1] Reject values outside optimization range
global rangemin
global rangemax
if min(x(:))<rangemin || max(x(:))>rangemax     %Outside  parameter range
    f=1e6*max([max(x) abs(min(x))]);
    %f=inf;           
    return
end


%[2] Unscaling CMAES parameters to  MODFLOW parameters
global OptRange
global ParaRange
xi=NaN(length(x),1);
for m=1:length(x)
    xi(m,1) = interp1(OptRange,ParaRange,x(m,1));
end

%[3.1] Transmisivity for Layer 1
HC1=reshape(xi(1:81),9,9);
fidT2=fopen('Tran2.dat','w');
for n=1:9
fprintf(fidT2,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n',HC1(n,:));
end
fclose(fidT2);

%[3.2] Transmisivity for Layer 2
HC2=reshape(xi(82:end),9,9);
fidT4=fopen('Tran4.dat','w');
for n=1:9
fprintf(fidT4,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n',HC2(n,:));
end
fclose(fidT4);

%[4] Run MODFLOW
!mf2005.exe twri.nam > NonConvergenceTest.txt
NC_Flag=NonConvergenceTest(); % NonConvergence test
if NC_Flag~=0                 % Return f=inf for  NonConvergence
    f=inf;
    return
end

%[5] Calculate objective function 
Hfid=fopen('head.hed','r');
SimuHead=fscanf(Hfid,'%f\n');
fclose(Hfid);
SimuHead=reshape(SimuHead,9,18)';
global TrueHead
f=sqrt(sum(sum((TrueHead-SimuHead).^2)));

