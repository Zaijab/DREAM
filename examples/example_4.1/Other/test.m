%0.013

tic
%[5] Calculate objective function 
Hfid=fopen('head.hed','r');
SimuHead=fscanf(Hfid,'%f\n');
fclose(Hfid);
SimuHead=reshape(SimuHead,9,18)';
load TrueHead
f=sqrt(sum(sum((TrueHead-SimuHead).^2)));
toc