function f=MODFLOW(HC,k)
close all
global t
global SMC
global MSMC
global MCMC
global fSMC
global fMSMC
global fMCMC
global MCMCr
%[3] Transmisivity for Layer 1
fidT2=fopen('Tran2.dat','w');
for n=1:9
    fprintf(fidT2,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n',HC(n,:));
end
fclose(fidT2);

%[4] Run MODFLOW
!mf2005.exe twri.nam > NonConvergenceTest.txt
NC_Flag=NonConvergenceTest(); % NonConvergence test
if NC_Flag~=0                 % Return f=inf for  NonConvergence
    f=inf;
    SMC(t,:,k) = SMC(t-1,:,k);
    MSMC(t,:,k) = MSMC(t-1,:,k);
    MCMC(t,:,k) = MCMC(t-1,:,k);
    MCMCr(t,:,k) = MCMCr(t-1,:,k);
    fSMC(t,k)=-inf;                  %An indication of error 2  
    fMSMC(t,k)=fMSMC(t-1,k);         %Previous f(x)
    fMCMC(t,k)=fMCMC(t-1,k);         %Previous f(x)
    return
end

%[5] Calculate objective function
Hfid=fopen('head.hed','r');
SimuHead=fscanf(Hfid,'%f\n');
fclose(Hfid);
SimuHead=reshape(SimuHead,9,18)';
load TrueHead
f=sqrt(sum(sum((TrueHead-SimuHead).^2)));



