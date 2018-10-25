function udot = mtlb_fRhs(t,u,s,data)

udot = zeros(5,1);

P = data.P;
Ep = data.Ep;
Imax  = data.Imax;
Sumax = data.Sumax;
Qsmax = data.Qsmax;
aE = data.aE;
aF = data.aF;
aS = data.aS;
Kf = data.Kf;
Ks = data.Ks;

Si = u(1);
Su = u(2);
Sf = u(3);
Ss = u(4);

Precip = P(s);
if Imax > 0
    EvapI  = Ep(s)*expFlux(Si/Imax,50);
    P_e    = P(s)*expFlux(Si/Imax,-50);
    Ep_e   = max(0.,Ep(s)-EvapI);
else
    EvapI = 0;
    P_e = P(s);
    Ep_e = Ep(s);
end
Evap   = Ep_e*expFlux(Su/Sumax,aE);
Perc   = Qsmax*expFlux(Su/Sumax,aS);
Runoff = P_e*expFlux(Su/Sumax,aF);
FastQ  = Sf/Kf;
SlowQ  = Ss/Ks;

udot(1) = Precip - EvapI - P_e;
udot(2) = P_e - Evap - Perc - Runoff;
udot(3) = Runoff - FastQ;
udot(4) = Perc - SlowQ;
udot(5) = FastQ + SlowQ;

end

function Qr = expFlux(Sr,a)
% Relative flux from exponential storage-flux relation
Sr = max(0.0,min(1.0,Sr));
if abs(a) < 1e-6
    Qr = Sr;    % approximately linear
else
    Qr = (1.-exponen(-a*Sr))/(1.-exponen(-a));
end
end

function f = exponen(x)
% Exponential function with protection against overflow
f = exp(min(300,x));
end
