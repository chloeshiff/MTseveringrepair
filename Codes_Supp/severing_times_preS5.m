%random walk model for severing competition
%stores many severing times for tubulin concentration 0,2,5 
%in files 'sevtimesct0.dat','sevtimesct2.dat','sevtimesct5.dat'
clear all;
close all;
MaxT=1000000;
numTraj=9000000;%make it high so we can reach 1000
%parameters
kt=.11; 
ks=.54; 
L=25;
cts=.01:.5:6.01;
for a=1:length(cts)
    ct=cts(a); %ct is concentration of tubulin
    sevtimes=zeros(1,1000);% to store distribution
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
     
        if sever==true
            sevdata=sevdata+1; %one additional severing event
            sevtimes(sevdata)=t(i);
        end
        if sevdata==1000 %we compute up to 1000 data points
            break
        end
    end %for trajectories

    %store the severing times as data files
    G=sevtimes;
    fileID4=fopen("sevtimesct"+num2str(ct)+".dat",'w');
    fprintf(fileID4,'%12.8f \r\n',G);
    fclose(fileID4);

end








