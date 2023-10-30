% S5, plots lvt, histogram of lengths, mean and cv of lengths for severing
% with repair, with repair and delay in severing time, and without repair
%requires severing times already computed (via severing_times_preS5) or the
%sevtimes ct files already in folder

histogramlvt=false; %true to plot histogram and lvt, set false to plot mean,cv
 
MaxT=1000000;
numTraj=100000;%for prob
numTraj=10;
ks=.54;
kt=.11;
x0=1;
L=25;
if histogramlvt
    cts=[0,2,5];
else
    cts=.01:1:5.01;
end
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
sevprobs=((r).^(x0)-1)./( ((r).^L)-1); %probability
mean1=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0));
means=zeros(1,length(cts));
stdevs=zeros(1,length(cts));
sst=1000; %time at which steady state is reached
stop_time=3000;% time simulation stops
cap=450000;
for n=[1,2,3]
    sevlens=zeros(2000,length(cts));
    ldiff=zeros(2000,length(cts));
    sevts=zeros(2000,length(cts));
    len_post_sev=zeros(2000,length(cts));
        if n==1
            comp=true;
            delay=true;
            titlen='with repair,with delay';
        elseif n==2
            comp=true;
            delay=false;
            titlen='with repair,w/o delay';
        elseif n==3
            comp=false;
            delay=false;
            titlen='w/o repair';
        end
    for a=1:length(cts)        
        sevdata=0;
        lengts=0;
        sevprob=sevprobs(a);
        if comp==false
           sevprob=1;
        end
        ct=cts(a);
        if histogramlvt
            filename="sevtimesct"+num2str(ct)+".dat";
            sevtimes=importdata(filename);
            ct=ct+1e-05;%need to add a bit to be nonzero but dont want to affect it
        end
        lengths=zeros(1,numTraj);
        c=1;
        sst_lengths=zeros(1,(stop_time-sst+1)*numTraj*20);
        num_proteins=zeros(1,(stop_time-sst+1)*numTraj*20);
        for j=1:numTraj
            l=zeros(1, MaxT);
            t=zeros(1, MaxT);
            dt=zeros(1, MaxT);
            numsev=0;
            severing=zeros(1,170000);
            for i=1:MaxT
                Nfree=ct;
                sev=0;
                k2=cs*(l(i)-1)*4/13*kds; %multiply by length of mt to multiply likelihood of landing by that of severing
                k1=ct*kdt;
                if l(i)<=0
                    k2=0; %cant sever if filament has 0 length
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
                    num_proteins(c)=nnz(severing);
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
                                if delay==false
                                    l(i+1)=s; %sever immediately
                                    sevdata=sevdata+1;
                                    lengts=lengts+(l(i)-s);
                                    sevts(sevdata,a)=t(i);
                                    sevlens(sevdata,a)=l(i);
                                    len_post_sev(sevdata,a)=l(i+1);
                                    break;
                                else
                                    sev=s;
                                    landlen=l(i);
                                    break;
                                end
                            end
                            c1=c1+w;
                            c2=c2+w;
                        end
                        sevtime=0;
                        if delay==true && severing(sev)==0 %no severing protein there
                             psevt=randi([1 1000]);
                             sevtime=sevtimes(psevt); 
                             severing(sev)=t(i+1)+sevtime;
                        end
                    end
                end
                severed=false;
                if delay
                    for f=1:l(i+1) %for any previous severings waiting to happen
                        if severing(f)>0 && t(i+1)>=severing(f) %if a severing is waiting and the time needed for severing has elapsed
                            lengt=f; 
                            l(i+1)=lengt;  %new length is where severing protein was
                            severing(f)=0; %severing no longer waiting at this position
                            %also get rid of all later severing proteins
                            if i>5000
                            lengts=lengts+(l(i)-lengt);
                            sevdata=sevdata+1;
                            sevlens(sevdata,a)=l(i);
                            ldiff(sevdata,a)=landlen-l(i);
                            len_post_sev(sevdata,a)=l(i+1);
                            sevts(sevdata,a)=t(i);
                            end
                            severed=true;
                            
                            break % lowest one is new length, don't need to go throught them all
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
    if histogramlvt
        figure(4)
        plot(t(1:i),l(1:i))
        hold on
        figure(3)
        histogram(sst_lengths(1:c-1))
        hold on
    else
        figure(1)
        plot(cts,means)
        figure(2)
        plot(cts,stdevs./means)
    end
end
figure(1)
% loglog(cts,th_means*4/13)
xlabel('[T] (\muM)');
ylabel('mean steady state length (nm)');
%legend;
saveas(gca,'meanssl.pdf')

figure(2)
xlabel('[T] (\muM)');
ylabel('CV (stdev/mean)');
%legend;
ylim([0 inf])
saveas(gca,'CVssl.pdf')

% if delay
    save('delayvarsbiggrowth','sevlens','sevts','cts','means','stdevs','len_post_sev','num_proteins');
