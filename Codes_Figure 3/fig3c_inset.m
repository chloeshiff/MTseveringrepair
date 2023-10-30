%plots theoretical mean-has the mean approxes that is in fig 3c inset
clear all;
close all
set(groot,'defaultLineLineWidth',1.5)
set(groot,'defaultAxesFontSize',13)
figure
MaxT=1000000;
cap=999090000; %parameter that specifies the max length we sum over in our probability calculation, must be very high for large ct, but simulation will take long
ks=.54;
kt=.11;
x0=1;
L=25;
cts=.01:.5:12.01; % vector of tubulin concentrations
th_means=zeros(size(cts));
thmean2=zeros(size(cts));
cs=.05;%.05
kds=12/6000; %
kdt=3.4;
pt=kt*cts./(kt*cts+ks);
ps=(ks./(kt*cts+ks));

r=pt./ps;
one_vec=1;

tstep=zeros(size(cts));
    %syms j
    for a=1:length(cts)
        tstep(a)=1/(ks+kt*cts(a));
    end
sevprobs= ((r).^(x0)-1)./( ((r).^L)-1); %probability
k=x0;
low_end=find(cts>6,1);
high_start=find(cts>4,1);
lows=zeros(1,low_end);
highs=zeros(1,length(cts)-high_start);
 for n=[2,3]
        if n==2
            comp=true; %with repair
            titlen='with repair';
        else
            comp=false; %without repair
            titlen='w/o repair';
        end
    th_means=zeros(size(cts));
    for a=1:length(cts)
        sevprob=sevprobs(a);
        if comp==false
           sevprob=1;
        end
        ct=cts(a);
        ra=ct*kdt;
        se=cs*kds*sevprob*4/13;
        rds=ra/se;
        th_dist=zeros(1,cap);
        %compute probability distribution values for the xes
        for leng=1:length(th_dist)
                th_dist(leng)=((leng)*se/ra)*exp(-(se/ra)*((leng)^2)/2);   
        end 
              th_dist=th_dist/sum(th_dist); %normalize
               th_dist_x=th_dist;
            %sum over probabilities times lengths to approximate mean
            %length
            for i=1:length(th_dist_x)
                th_means(a)=th_means(a)+(i-1)*th_dist_x(i);    
            end
        thmean2(a)=sqrt(pi*rds/2);
    %compute low [T] and high [T] approximations
        if a<=low_end
            if comp
                lows(a)=sqrt(pi*kdt*ct/(2*kds*cs*4/13*(1-kt*ct/ks)));
            else
                lows(a)=sqrt(pi*kdt*ct/(2*kds*cs*4/13));
            end
        end
        if a>=high_start
            if comp
            highs(a-high_start+1)=sqrt(pi*kdt*ct/(2*kds*cs*4/13*(ks/(kt*ct))^(L-1)));
            else
                highs(a-high_start+1)=sqrt(pi*kdt*ct/(2*kds*cs*4/13));
            end
        end
    end
    %plot mean on loglog plot
    loglog(cts,th_means*4/13/1000)
    hold on
        %plot approximations for severing with repair
        if n==2
            loglog(cts(1:low_end),lows*4/13/1000, '-','Color','black');
            hold on
            loglog(cts(high_start:length(cts)),highs*4/13/1000, '-','Color','black');
            hold on
        end
 end
    
   