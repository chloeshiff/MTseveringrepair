%Figure S2
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',15)
figure;
cts=[0,2,5];


% plot points
data=readmatrix("vemu_dat_fromantonina.xlsx");
T0=data(:,1);
T0=T0(~isnan(T0));
T2=data(:,2);
T2=T2(~isnan(T2));
T5=data(:,3);
T5=T5(~isnan(T5));

%histograms of data

figure
histogram(T0,NumBins=10,Normalization="probability");
simct0=load('sevtimesct0.dat');
hold on
histogram(simct0,Normalization="probability")


figure
histogram(T2,NumBins=10,Normalization="probability");
simct2=load('sevtimesct2.dat');
hold on
histogram(simct2,Normalization="probability")


figure
histogram(T5,NumBins=9,Normalization="probability");
simct5=load('sevtimesct5.dat');
hold on
histogram(simct5,Normalization="probability")
