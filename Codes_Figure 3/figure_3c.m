% plots main figure and inset for figure 3c, which is mean and cv
MaxT=1000000;
numTraj=140;
ks=.54;
kt=.11;
x0=1;
L=25;
cts=.01:5.51;
th_means=zeros(1,length(cts));
cs=.05;
kds=12/6000;
kdt=3.4;
pt=kt*cts./(kt*cts+ks);
ps=(ks./(kt*cts+ks));
tstep=zeros(size(cts));
for a=1:length(cts)
    tstep(a)=1/(ks+kt*cts(a));
end
r=pt./ps;
sevprobs=((pt./ps).^(x0)-1)./( ((pt./ps).^L)-1); %probability
mean1=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0));
means=zeros(1,length(cts));
stdevs=zeros(1,length(cts));
sst=9000; %time at which steady state is reached
stop_time=12000;
cap=450000;
for n=[2,3]
        if n==2
            comp=true;
            titlen='with repair';
        else
            comp=false;
            titlen='w/o repair';
        end
    for a=1:length(cts)
        sevprob=sevprobs(a);
        if comp==false
           sevprob=1;
        end
        ct=cts(a);
        lengths=zeros(1,numTraj);
        c=1;
        sst_lengths=zeros(1,(stop_time-sst+1)*numTraj*20);
        for j=1:numTraj
            l=zeros(1, MaxT);
            t=zeros(1, MaxT);
            dt=zeros(1, MaxT);
            numsev=0;
            severing=zeros(1,70000);
            for i=1:MaxT
                seva=0; %reset amount severed to 0 for each step
                Nfree=ct;
                sev=0;
                k2=cs*(l(i)-1)*4/13*kds; %multiply by length of mt to multiply likelihood of landing by that of severing
                k1=ct*kdt;
                if l(i)<=0
                    k2=0;
                end
                k0=k1+k2;
                r1=rand;
                dt(i)=(1/k0)*log(1/r1);
                if t(i)>stop_time
                    lengths(j)=l(i);
                    break;
                end
                if t(i)>sst
                    sst_lengths(c)=l(i);
                    c=c+1;
                end
                t(i+1)=t(i)+dt(i);
                p=rand;
                if p<=(k1/k0)
                    l(i+1)=l(i)+1;
                else
                    l(i+1)=l(i); %length remains until a change occurs- in this time step only a severing protein lands
                    sevrand=rand;
                    if sevrand<sevprob
                        w=(1-(k1/k0))/(l(i)-1);  %break up remaining possible p into segments
                        c1=(k1/k0); %lowest possible random variable
                        c2=c1+w; %increment for each possible segment
                        for len=1:l(i) %for every possible length
                            s=len;
                            if p>c1 && p<=c2 %if p is in corresponding interval, land at this position
                                if comp==false
                                    l(i+1)=s; %no longer sever immediately
                                    break;
                                else
                                    sev=s;
                                    break;
                                end
                            end
                            c1=c1+w;
                            c2=c2+w;
                         end
                        if comp==true && severing(sev)==0 %no severing protein there
                            sevtime=mean1(a);
                            severing(sev)=t(i+1);
                        end
                    end
                end
                severed=false;
                if comp
                    for f=1:l(i+1) %for any previous severings waiting to happen
                        if severing(f)>0 && t(i+1)>=severing(f) %if a severing is waiting and the time needed for severing has elapsed
                            lengt=f; %want to sever maximum amount on each step
                            l(i+1)=lengt;  %sever this amount from current length, not allowing less than 0
                            severing(f)=0; %severing no longer waiting at this position
                            %also get rid of all later severing proteins
                            severed=true;
                            break
                        end
                    end
                    if severed==true
                        lastsev=find(severing,1,'last'); %gives index of last(highest position) severing waiting
                        for m=lengt+1:lastsev
                            severing(m)=0;
                        end
                    end
                end
            end
        end
        means(a)=mean(sst_lengths(1:c-1));
        stdevs(a)=std(sst_lengths(1:c-1));
        % compute theory
        rds=(ct*kdt)/(cs*kds*sevprob);
        th_means(a)=exp(rds)*rds^(1-rds)*(gamma(rds)-gammainc(rds,rds));
        ct=cts(a);
        ra=ct*kdt;
        se=cs*kds*sevprob*4/13;
        rds=ra/se;
        th_dist=zeros(1,cap);
        for leng=1:length(th_dist)
                th_dist(leng)=(leng*se/ra)*exp(-(se/ra)*(leng^2)/2);   
        end 
            th_dist_x=th_dist*4/13/1000;
            for i=1:length(th_dist_x)
                th_means(a)=th_means(a)+i*th_dist_x(i);             
            end
    end
 
    figure(1)
    loglog(cts,means*4/13/1000,'o','DisplayName',titlen)
    hold on
    loglog(cts,th_means)
    hold on
    disp(means)
    figure(2)
    plot(cts,stdevs./means,'DisplayName',titlen)
    disp(sum(stdevs./means)/length(means))
    hold on

end
figure(1)
xlabel('[T] (\muM)');
ylabel('mean steady state length (nm)');

figure(2)
xlabel('[T] (\muM)');
ylabel('CV (stdev/mean)');
%legend;
ylim([0 inf])



