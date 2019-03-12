clear all
clc
load OD162_L50
[dummy indx]=min(outdat(:,1));
for gen=0:100:indx
    if gen==0
        continue
    end
    if outdat(gen,1)<200
        x=outdat(gen,3:end)';
        PlotGen(x,gen,outdat(gen,1));
    end
end
