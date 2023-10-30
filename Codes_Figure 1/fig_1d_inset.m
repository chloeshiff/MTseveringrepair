%% setup
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',15)
load('sevtimestderr_for_fig1.mat') %from standard_deviation_sevtime_fig1.m
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
    bootstdev=bootstrp(numTraj2,@std,T); %draws 200 sample severing times, computes stdev
    means(a)=mean(bootmean);
    stdevs(a)=std(bootmean);
    mstd(a)=mean(bootstdev);
    mstderr(a)=mean(bootstdev)/sqrt(numTraj2);
    stdstd(a)=std(bootstdev);
    stderr(a)=std(bootstdev)/sqrt(numTraj2);
end

figure
errorbar(cts,mstd,stdstd,'o','LineWidth',1,'Color','blue')
hold on
%plot(cts,data(:,2),'-o')
hold on
%RUN A CODE FIRST??- standard_deviation_sevtime.m nope that code saves and
%this loads things 
cts2=0:.5:5;
errorbar(cts2,sevtimesstdmean,sevtimesstd2,'-o','DisplayName','simulation','Color','red','LineWidth',1);
xlabel('[T]')
saveas(gca,'std_bootstrapping.pdf')