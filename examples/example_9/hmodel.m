function SimQ = hmodel(Pars,Extra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function runs conceptual rainfall-runoff model (hmodel)
%INPUT
%Pars: column vector of model parameters
%Extra: structure containing boundary conditions
%
%OUTPUT
%SimQ: column vector of simulated discharge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract model parameters
par = Extra.fpar;            %fixed parameters
par(Extra.idx_vpar) = Pars;  %variable parameters
par = par';                  %make it a column vector
detpar = par(1:end-11);

%Specify dimensions and settings for numerical integration
nsubc = 1;
npar = length(detpar);
ns = length(Extra.Precip);
ds = 1;
dtmin = 1e-0;
dtmax = 1e-0;
dtini = dtmax;
abstol = 1e-3;
reltol = 0;
imethod = 6;
mo = [1 1 1 4 1 1];
m_order = mo(imethod);
Specs = [nsubc npar ns ds dtmin dtmax dtini abstol reltol imethod m_order]';

%Initial states
States_in = [1e-10 1e-10 1e-10 1e-10]';

%Run model
A_in = [Specs; Extra.Precip; Extra.Ep; detpar; States_in];
A_out = hmodel_mat(A_in);

%Get output
SimQ = A_out(1:ns);
SimQ = SimQ(Extra.idx);

end