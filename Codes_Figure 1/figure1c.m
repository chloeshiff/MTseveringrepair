%%
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',15)
%parameters
kt=.11; 
ks=.54;  % rate per second of a tubulin protein falling off
L=25;
cts=0:.01:5;
x0=1;
%%
kT=kt*cts;
pt= kT./(kT+ks);
ps= ks./(kT+ks);
r=pt./ps;
tstep=zeros(size(cts));
one_vec=ones(size(cts));
for a=1:length(kT)
    tstep(a)=1/(ks+kT(a));
end
C5=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0)); %severing time
for i=1:length(kT)
    if kT(i)==ks
        C5(i)=(L^2-x0^2)/3;
    end
end
figure
plot(cts, C5,'DisplayName','theory')
%% setup
load('sevtimestderr_for_fig1.mat')
data=readmatrix("vemu_dat_fromantonina.xlsx");
T0=data(:,1);
T0=T0(~isnan(T0));
T2=data(:,2);
T2=T2(~isnan(T2));
T5=data(:,3);
T5=T5(~isnan(T5));
cts=[0,2,5];
means=zeros(1,length(cts));
stdevs=zeros(1,length(cts));
mstd=zeros(1,length(cts));
mstderr=zeros(1,length(cts));
stdstd=zeros(1,length(cts));
stderr=zeros(1,length(cts));
data=[32.3,9.2; %0: mean, std 
      75.4,30.8; %2: mean, std 
      193.8,101.5] ;

%%
for a=1:length(cts)
    if a==1
        T=T0;
    elseif a==2
        T=T2;
    else
        T=T5;
    end
    
    bootmean= bootstrp(numTraj2,@mean,T);
    bootstdev=bootstrp(numTraj2,@std,T); %draws sample severing times, computes stdev
    means(a)=mean(bootmean);
    stdevs(a)=std(bootmean);
    mstd(a)=mean(bootstdev);
    mstderr(a)=mean(bootstdev)/sqrt(numTraj2);
end
hold on
errorbar(cts,data(:,1),mstderr,'o','DisplayName','experiment w bootstrapping')
xlabel('[T]')
ylabel('mean severing time')
legend
