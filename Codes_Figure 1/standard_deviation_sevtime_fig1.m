%random walk model for severing competition
%plots mean end time and severing time vs concentration of tubulin
set(groot,'defaultLineLineWidth',1.8)
set(groot,'defaultAxesFontSize',15)
clear all;
%close all;
figure;
MaxT=1000000;
numTraj2=750;
numTraj=20000;

%parameters
kt=.11; 
ks=.54;  % probability(ish-actually a rate) of a tubulin protein falling off
L=25;
cts=0:.5:5;
numdata=[63,34,27];
probs=zeros(numTraj2,length(cts));
sevtimes=zeros(numTraj2,length(cts));
sevtimesstd=zeros(numTraj2,length(cts));
sevtimesstderr=zeros(numTraj2,length(cts));
endtimes=zeros(numTraj2,length(cts));

meanendtime=zeros(1,length(cts));
stdevendtime=zeros(1,length(cts));
meansevtime=zeros(1,length(cts));
stdevsevtime=zeros(1,length(cts));

tstep=zeros(1,length(cts));
prob=zeros(1,length(cts));
 psevering=zeros(1,length(cts));
 x0=1;
for m=1:numTraj2
    for a=1:length(cts)
        ct=cts(a); %ct is concentration of tubulin
        endtime=zeros(1,numTraj);
        sevtime=zeros(1,numTraj);
        ttosever=0;
        sevdata=0; %to store number of severing events
        for j=1:numTraj
            t=zeros(1,MaxT);
            dt=zeros(1,MaxT);
            l=zeros(1,MaxT);
            l(1)=x0;
            mt=x0;
            for i= 1:MaxT
                kp=ct*kt; % prob of tubulin adding, based on concentration
                k0=ks+kp;%rate of something happening is sum of rate of each indepedent event occurring
                r1=rand;
                dt(i)=(1/k0)*log(1/r1);
                 if mt==0  %if all tubulin spots filled 
                     sever=false;
                     break;
                 end
                if mt==L %if all tubulin gone
                    sever=true;
                     break;
                end
                t(i+1)=t(i)+dt(i);
                p=rand;
                if p<=(ks/k0)  %a tubulin comes off- go towards severing
                  mt=mt+1;
                  l(i+1)=l(i)+1;
                else %repair-fill one unit
                    mt=mt-1;
                    l(i+1)=l(i)-1;
                end
            
            end %for i
            endtime(j)=t(i);
            if sever==true
                sevdata=sevdata+1; %one additional severing event
                ttosever=ttosever+t(i); %add to total time
                sevtime(sevdata)=t(i); %instead add to a vector to make stdec easier
            end
%             if sevdata==numdata(a) %sample as many times as in vemu
%                     break;
%             end
                if sevdata>1000
                    break;
                end
        end %for trajectories
        sevtime=sevtime(1,1:sevdata); %cut off values that might be 0 at end so we just have a vector of all the times it took to sever
        prob(a)=sevdata/numTraj; %fraction successful severing events
        meanendtime(a)=mean(endtime);%store values for each ct
        stdevendtime(a)=std(endtime);
        %meansevtime(a)=ttosever/sevdata; %sum of severing times/num severing events
        meansevtime(a)=mean(sevtime);
        stdevsevtime(a)=std(sevtime);

        psevtimes=sevtimes./sum(sevtimes); %turns frequency into probability
       
    end
   probs(m,:)=prob; %store values for each ct for each trial
   sevtimes(m,:)=meansevtime;
   sevtimesstd(m,:)=stdevsevtime; %store std for each trial
   sevtimesstderr(m,:)=stdevsevtime/sqrt(sevdata); %store std for each trial
   endtimes(m,:)=meanendtime;
end
prob2=mean(probs);
probstd=std(probs);
sevtime2=mean(sevtimes);%mean of mean severing times
sevtimestd=std(sevtimes); %stdev of mean severing times
endtime2=mean(endtimes);
endtimestd=std(endtimes);
sevtimesstdmean=mean(sevtimesstd); %mean of stds acroaa trials
sevtimesstd2=std(sevtimesstd); %take std of the stds across trials
sevtimesstderr2=std(sevtimesstderr)/sqrt(numTraj2);
%%
kT=kt*cts;
pt= kT./(kT+ks);
ps= ks./(kT+ks);
r=pt./ps;
k=ones(size(cts)).*x0;
n=ones(size(cts)).*L;
one_vec=ones(size(cts));
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
C5=tstep.*(((r+one_vec)./(r-one_vec)).*(((r.^n+one_vec)./(r.^n-one_vec)).*n-((r.^k+one_vec)./(r.^k-one_vec)).*k));
            for i=1:length(kT)
                if kT(i)==ks
                    C5(i)=(L^2-x0^2)/3;
                end
            end

%%
figure;
%errorbar(cts,prob2,probstd,'-o','DisplayName','simulation');
%xlabel('concentration of tubulin');
%ylabel('probability of severing');
%title("kt="+num2str(kt)+"  "+"ks="+num2str(ks)+" "+"x0="+num2str(x0)+" "+"L="+num2str(L));
hold on;
plot(cts,C,'DisplayName','theory');
hold on;
%plot(cts,B4,'DisplayName','theory3');
legend;
saveas(gca,'sevproberr.pdf');

figure;
errorbar(cts,sevtimesstdmean,sevtimesstd2,'-o','DisplayName','simulation');
% xlabel('concentration of tubulin');
% ylabel('mean severing time');
% title("kt="+num2str(kt)+"  "+"ks="+num2str(ks)+" "+"x0="+num2str(x0)+" "+"L="+num2str(L));
hold on;
cts2=[0,2,5];
plot(cts2,data(:,2),'o','DisplayName','experiment');
legend;
saveas(gca,'sevtimestderr.pdf');

save('sevtimestderr_for_fig1.mat','sevtimesstdmean','sevtimesstd2','sevtimesstderr2','numTraj2');

%%%sevtime for many trials, stdev for each trial, mean of stds across the trials for std,
%%%std across the trials over sqrt numtrials gives std?