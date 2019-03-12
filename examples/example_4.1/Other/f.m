function fitness=f(xStar,k)

global x
global t
global f
global fBest
global fx

%[2] Metropolis sampler
fStar=f(xStar);                     %New x
fPrev=fx(t-1,k);                    %Previous x
alpha = min([1, (fStar/fPrev)*1]);  %Acceptance ratio

if rand < alpha             %Accept
    x(t,:,k) = xStar;
    fx(t,k)=fStar;
else                        %Reject
    x(t,:,k) = x(t-1,:,k);
    fx(t,k)=fPrev;
end

%[3] Fitness
fitness=-alpha;
if fStar>fBest(1,k)
    fBest(1,k)=fStar;
end

