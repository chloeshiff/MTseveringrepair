%random walk model for severing competition

%parameters
kt=.11; 
ks=.54;  % rate per second of a tubulin protein falling off
L=25;
cts=0:.01:30;
x0=1;
%%
kT=kt*cts;
pt= kT./(kT+ks);
ps= ks./(kT+ks);
r=pt./ps;
tstep=zeros(size(cts));
for a=1:length(kT)
    tstep(a)=1/(ks+kT(a));
end
C= ((pt./ps).^(x0)-1)./( ((pt./ps).^L)-1); %probability
data=[32.3,9.2; %0: mean, std 
      75.4,30.8; %2: mean, std 
      193.8,101.5] ;
C5=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0)); %severing time
for i=1:length(kT)
    if kT(i)==ks
        C5(i)=(L^2-x0^2)/3;
    end
end
%% loglog plot for probability
figure
loglog(cts,C,'DisplayName','theory');
hold on;
ct4=0:.005:4.895;
plinet=zeros(1,length(ct4));
for i=1:length(ct4)
    plinet(i)=1-(kt*ct4(i)/ks); %asymptotics
end
loglog(ct4,plinet,'--','Color','black')
ct5=5:.01:30;
pcurvet=zeros(1,length(ct5));
for i=1:length(ct5)
    pcurvet(i)=(kt*ct5(i)/ks)^(-24); %our asymptotic approx
end
loglog(ct5,pcurvet,'--','Color','black')
%% plot for severing time
figure;
hold on;
plot(cts,C5,'DisplayName','theory');
hold on
ct2=0:.01:4;
linet=zeros(1,length(ct2));
for i=1:length(ct2)
    linet(i)=(L/ks)*(1+(kt*ct2(i)/ks));
end
plot(ct2,linet,'--','Color','black')
hold on
ct3=12:.01:30;
curvet=zeros(1,length(ct3));
for i=1:length(ct3)
    curvet(i)=(L/ks)*(1/(kt*ct3(i)/ks)); %our asymptotic approx
end
plot(ct3,curvet,'--','Color','black')
%% plots [T] vs probability
figure
plot(cts,C)
%% loglog plot for severing time
figure;
loglog(cts,C5,'DisplayName','theory');
hold on
ct2=0:.01:4;
linet=zeros(1,length(ct2));
for i=1:length(ct2)
    linet(i)=(L/ks)*(1+(kt*ct2(i)/ks)); %asymptotic
end
loglog(ct2,linet,'--','Color','black')
hold on
ct3=12:.01:30;
curvet=zeros(1,length(ct3));
for i=1:length(ct3)
    curvet(i)=(L/ks)*(1/(kt*ct3(i)/ks));
end
loglog(ct3,curvet,'--','Color','black')