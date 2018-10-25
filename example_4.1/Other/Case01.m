function f=Case01(varargin)

global CMARange
global ParaRange
global rangemin
global rangemax
global t
global SMC
global MSMC
global MCMC
global fSMC
global fMSMC
global fMCMC
global MCMCr
fclose('all');

%[1] Reject values outside optimization range
x=varargin{1,1};
if min(x(:))<rangemin || max(x(:))>rangemax     %Outside  parameter range
    f=1e6*max([max(x) abs(min(x))]);  
    k=varargin{1,2};
    SMC(t,:,k) = SMC(t-1,:,k);
    MSMC(t,:,k) = MSMC(t-1,:,k);
    MCMC(t,:,k) = MCMC(t-1,:,k);
    MCMCr(t,:,k) = MCMCr(t-1,:,k);
    fSMC(t,k)=inf;               %An indication of error 1           
    fMSMC(t,k)=fMSMC(t-1,k);     %Previous f(x)
    fMCMC(t,k)=fMCMC(t-1,k);     %Previous f(x)
    return
end

%[2] Unscaling CMAES parameters to  MODFLOW parameters
xi= interp1(CMARange,ParaRange,x);
HC=reshape(xi,9,9);
k=varargin{1,2};
f=MODFLOW(HC,k);

%[2]Samplers
%(2.1) SMC sampler
fStar=f;
SMC(t,:,k) = xi';
fSMC(t,k)=fStar;

%(2.2) Metropolis SMC sampler
fPrev=fMSMC(t-1,k);                 %Previous x
alpha = min([1, (fStar/fPrev)*1]);  %Acceptance ratio
if rand < alpha             %Accept
    MSMC(t,:,k) = xi';
    fMSMC(t,k)=fStar;
else                        %Reject
    MSMC(t,:,k) = MSMC(t-1,:,k);  
    fMSMC(t,k)=fPrev;
end

%(2.3) Metropolis MCMC sampler
e=varargin{1,3};                            %CMA Range
xr=MCMCr(t-1,:,k)+e';                       %CMA Range
xStar= interp1(CMARange,ParaRange,xr);      %MF Range
HC=reshape(xStar,9,9);
fStar=MODFLOW(HC,k);
fPrev=fMCMC(t-1,k);                 %Previous x
alpha = min([1, (fStar/fPrev)*1]);  %Acceptance ratio

if rand < alpha             %Accept
    MCMC(t,:,k) = xStar';
    MCMCr(t,:,k)=xr';
    fMCMC(t,k)=fStar;
else                        %Reject
    MCMC(t,:,k) = MCMC(t-1,:,k);  
    MCMCr(t,:,k)=MCMCr(t-1,:,k);
    fMCMC(t,k)=fPrev;
end



