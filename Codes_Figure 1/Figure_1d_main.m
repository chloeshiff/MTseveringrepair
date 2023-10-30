%FIGURE 1 CODE
%
%% 1D
data=readtable('vemu_dat_fromantonina.xlsx');
T0=data.Var1(~isnan(data.Var1));
T2=data.Var2(~isnan(data.Var2));
T5=data.Var3(~isnan(data.Var3));
MaxT=1000000;
numTraj=400000;
kt=.11; 
ks=.54;  % probability(ish-actually a rate) of a tubulin protein falling off
L=25;
cts=[0,2,5];
meanendtime=zeros(1,length(cts));
stdevendtime=zeros(1,length(cts));
meansevtime=zeros(1,length(cts));
stdevsevtime=zeros(1,length(cts));
tstep=zeros(1,length(cts));
prob=zeros(1,length(cts));

%% simulation and plot histogram
for a=1:length(cts)
    ct=cts(a); %ct is concentration of tubulin
    endtime=zeros(1,numTraj);
    sevtime=zeros(1,numTraj);
    leftsev=zeros(1,numTraj);
    ttosever=0;
    sevdata=0; %to store number of severing events
    for j=1:numTraj
        t=zeros(1,MaxT);
        dt=zeros(1,MaxT);
        l=zeros(1,MaxT);
        left=0;
        x0=1;
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
                left=left+1;
            end

        end %for i
        endtime(j)=t(i);
        if sever==true
            sevdata=sevdata+1; %one additional severing event
            ttosever=ttosever+t(i); %add to total time
            sevtime(sevdata)=t(i); %instead add to a vector to make stdec easier
        end
     end %for trajectories
    sevtime=sevtime(1,1:sevdata); %cut off values that might be 0 at end so we just have a vector of all the times it took to sever
    figure;
    %hold on
    bins=[8,10,20];
    histogram(sevtime,'Normalization','Probability')
    hold on
    [N,edges]=histcounts(sevtime);
    binwidth=edges(2)-edges(1);
    xes=edges(2:end)-binwidth/2; %get middle of each bin
    plot(xes,N/length(sevtime),'black');
    fname="sevtimedistn"+num2str(ct)+"500.pdf";
    xlim([0 500])
    ylim([0 .12])
    saveas(gca,fname);
    set(gca,'color','none')
    xlim([0 500]);
    %xlabel('severing time (s)');
    ylim([0 .18]);
    %title("[T]="+num2str(ct)+" "+"\muM");
    prob(a)=sevdata/numTraj; %fraction successful severing events
    meanendtime(a)=mean(endtime);
    stdevendtime(a)=std(endtime);
    meansevtime(a)=mean(sevtime);
end
