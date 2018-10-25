function SimuHead =Case01(Pars,~)
global DREAM_dir EXAMPLE_dir
fclose('all');

%[1]Write input(HC for Layer 1)
cd(EXAMPLE_dir)
HC=reshape(Pars,9,9);
fidT2=fopen('Tran2.dat','w');
for n=1:9
    fprintf(fidT2,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n',HC(n,:));
end
fclose(fidT2);

%[2] Run MODFLOW
!mf2005.exe twri.nam > NonConvergenceTest.txt
NC_Flag=NonConvergenceTest(); % NonConvergence test
if NC_Flag~=0                 % Return Huge head values for  NonConvergence
    SimuHead=ones(9*9*2,1)*1000;
    cd(DREAM_dir)
    return
end

%[3] Read output (head)
Hfid=fopen('head.hed','r');
SimuHead=fscanf(Hfid,'%f\n');
SimuHead=reshape(SimuHead,9,18)';
SimuHead=SimuHead(:);
fclose(Hfid);
cd(DREAM_dir)








