% severing time theory
x0=1;
kt=.11; 
ks=.54;  % probability(ish-actually a rate) of a tubulin protein falling off
L=25;
cts=0:.2:10;
kT=kt*cts;
pt= kT./(kT+ks);
ps= ks./(kT+ks);
r=pt./ps;
for a=1:length(kT)
    tstep(a)=1/(ks+kT(a));
end

B= (x0./(1- 2*ps)) - (L./(1- 2*ps)).*(((pt./ps).^(x0)-1)./( ((pt./ps).^L)-1)); %normal boundary conditions
B=tstep.*B;
C= ((pt./ps).^(x0)-1)./( ((pt./ps).^L)-1); %probability
M=B./C;
data=[32.3,9.2; %0: mean, std 
      75.4,30.8; %2: mean, std 
      193.8,101.5] ;
C5=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0));
for i=1:length(kT)
    if kT(i)==ks
        C5(i)=(L^2-x0^2)/3;
    end
end

saveas(gca,'response_curve.pdf')

%%
%% plot mean severing time vs [T] for multiple Ls
figure;
Ts=0:.01:10;
x0=1;
means=[32.3,75.4,193.8];
Ls=[25,30,40];
kss=[.54,.7,.98];
kts=[.11,.15,.22];
for i=[1,2,3]
    L=Ls(i);
    ks=kss(i);
    kt=kts(i);
    kT=kt*cts;
    r=(kt/ks)*Ts;
    C52=(1./(kt.*Ts+ks*1)).*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0));
    C= ((r).^(x0)-1)./( ((r).^L)-1); %probability
    figure(1)
    plot(Ts,C52,'DisplayName',"N="+num2str(L));
    hold on
    figure(2)
    plot(Ts,C,'DisplayName',"N="+num2str(L));
    hold on
end
figure(1)
legend
figure(2)
legend
saveas(gca,'response_curve_multipleLs.pdf')

%%
x0s=[1,2,3];
Ls=[25,30,40];
kss=[.54,.7,.98];
kts=[.11,.15,.22];
for i=[1,2,3]
    L=Ls(i);
    x0=x0s(i);
    ks=kss(1);
    kt=kts(1);
    kT=kt*cts;
    r=(kt/ks)*Ts;
    C52=(1./(kt.*Ts+ks*1)).*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0));
    C= ((r).^(x0)-1)./( ((r).^L)-1); %probability
    figure(4)
    plot(Ts,C52,'DisplayName',"x_0="+num2str(x0));
    hold on
    figure(3)
    plot(Ts,C,'DisplayName',"x_0="+num2str(x0));
    hold on
end
figure(4)
legend
figure(3)
legend
saveas(gca,'response_curve_multipleLs.pdf')