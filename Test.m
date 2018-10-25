clc
clear all
close all
load Out.mat

%[4] Output
%(4.0) General opitions
nBins =[20 20];
ra=[-5 5];
rb=[-5 5];
a = linspace(ra(1), ra(2), nBins(1));
b = linspace(rb(1), rb(2), nBins(2));
[A B]=meshgrid(a,b);

%(4.1) Plot analytical density with proposal distribution
subplot(121);
f=@(x) mvnpdf(x,[0 0],[2 .9;.9 2]);  %Bivariate normal
PD = f([A(:), B(:)]);      %Probability density of x f(x,mu,Cov)
PD=reshape(PD,nBins(2),nBins(1));
surf(A,B,PD,'EdgeColor','none','FaceLighting','phong'); 
xlabel('A'); ylabel('B'); zlabel('p({\bfx})');
title('Analytic Distribution')
axis tight

%(4.2) Plot sampled distribution
subplot(122);
x=ParSet(:,1:2);
nS=size(x,1);
F=hist3(x, {a b});          %Frequency 
SF=F/trapz(b,trapz(a,F));   %Probability denisty of x
surf(A,B,SF','EdgeColor','none','FaceLighting','phong');
xlabel('A'); ylabel('B'); zlabel('p({\bfx})');
title('Sampled Distribution');
axis tight

%(4.3) Error
Err=mean(mean(abs(PD'-SF)));
disp(['Number of samples= ' num2str(nS) ' Error= ' num2str(Err/1*100) ' %'])




