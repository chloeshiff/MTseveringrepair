% bind time and severing time vs [T]
%severing time we already have
%bind time is a function of [T] and mean length
close all;
set(groot,'defaultLineLineWidth',1.8)
set(groot,'defaultAxesFontSize',15)
set(gca,'FontWeight','bold')
kdt=3.4;
cs=.05;
kds=12/6000;
cts=0:.1:5.5;
th_means=zeros(size(cts));
sevtimes=zeros(size(cts));
C5=zeros(size(cts));
bind_time=zeros(size(cts));
kt=.11;
ks=.54;
L=25;
x0=1;
for a=1:length(cts)
    ct=cts(a);
    sevprob=((kt*ct/ks)^(x0)-1)/(((kt*ct/ks)^L)-1); 
    ra=ct*kdt;
    se=cs*kds*sevprob*4/13; %for length with repair- a bit off since length differs with delay
    cap=20000;
    th_dist=zeros(1,cap);
    for leng=1:length(th_dist)
            th_dist(leng)=(leng*se/ra)*exp(-(se/ra)*(leng^2)/2);   
    end 
        th_dist_x=th_dist*4/13; % gives length in nm
        for i=1:length(th_dist_x)
            th_means(a)=th_means(a)+i*th_dist_x(i);   
        end
    tstep=kt*ct+ks;
    r=kt*ct/ks;
    C5(a)=tstep*(((r+1)/(r-1))*(((r^L+1)./(r^L-1))*L-((r^x0+1)/(r^x0-1))*x0));
            if kt*ct==ks
                C5(i)=tstep*(L^2-x0^2)/3;
            end
    bind_time(a)=1/(th_means(a)*(kds*cs)*sevprob);
end
plot(cts,bind_time,'DisplayName','bind time')
hold on
plot(cts, C5,'DisplayName','severing time')
legend;
xlabel('concentration of tubulin (\muM)')
ylabel('time (s)')