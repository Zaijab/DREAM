function [log_likelihood]=DCT_GPR(Pars, Extra)

for k=1:Extra.parz
    for i=1:Extra.parx
        Extra.model_DCT(k,i)=Pars((k-1)*Extra.parx+i); % Assign proposed model
    end
end                     
model=idct2(Extra.model_DCT); % Do inverse DCT
model=10.^model; %Estimated parameters are given in logarithmic scale;
model=model';
minvalue=min(min(model));
maxvalue=max(max(model));

data=Extra.J*model(:)-Extra.datasim; % Calculate residual
data=data./Extra.error;              % Calculate weighted residual
log_likelihood=sum(data.^2); % Calculate log-likelihood
if minvalue<(1/0.17)
   		log_likelihood=log_likelihood*(1+(1/0.17)/minvalue); % Penalty of models outside range
end

if maxvalue>(1/0.05)
        log_likelihood=log_likelihood*(1+maxvalue/(1/0.05)); % Penalty of models outside range
end
log_likelihood=-1/2*log_likelihood;