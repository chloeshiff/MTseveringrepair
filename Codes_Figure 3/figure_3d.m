% plots main figure for figure 3e, the length distributions
clear all
figure;
MaxT=1000000;
%numTraj=155;
numTraj=19;
ks=.54;
kt=.11;
x0=1;
L=25;
cts=[5];
th_means=zeros(size(cts));
cs=.05;
kds=12/6000;
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
sevprob= ((pt./ps).^(x0)-1)./( ((pt./ps).^L)-1); %probability
means1=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-one_vec)).*L-((r.^x0+one_vec)./(r.^x0-one_vec)).*x0));
stdev1=1;

means=zeros(1,length(cts));
stdevs=zeros(1,length(cts));
sevp=0; %2 if keep severing prob
sst=2500; %time at which steady state is reached
stop_time=6000;
 for n=[2,3]
        if n==1
            comp=true;
            delay=true;
            titlen='with repair,with delay';
        elseif n==2
            comp=true;
            delay=false;
            titlen='with repair,w/o delay';
        else
            comp=false;
            delay=false;
            titlen='w/o repair';
        end
    for a=1:length(cts)
        mean1=means1(a);
        ct=cts(a);
        prob=zeros(1,300000); %length is largest possible mt length (prediction)
        tavg=zeros(1,MaxT);
        lavg=zeros(1,MaxT);
        ivals=zeros(1,numTraj);
        lengths=zeros(1,numTraj);
        c=1;
        sst_lengths=zeros(1,(stop_time-sst+1)*numTraj*20);
        for j=1:numTraj
            l=zeros(1, MaxT);
            t=zeros(1, MaxT);
            dt=zeros(1, MaxT);
            numsev=0;
            %to store severings waiting- index is amount severed/position on filament
            % from left end and value is time needed to reach for severing to occur
            severing=zeros(1,length(prob));
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
                    if comp==false
                        sevprob=1;
                    end
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
                            sevtime=normrnd(mean1,stdev1);
                            if delay==false
                                sevtime=0;
                            end
                            severing(sev)=t(i+1)+sevtime;
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
            ivals(j)=i;
            for o=1:i
                tavg(o)=tavg(o)+t(o);
                lavg(o)=lavg(o)+l(o);
            end
        end
        mini=min(ivals);
        t2=zeros(1,mini);
        l2=zeros(1,mini);
        for h=1:mini
            t2(h)=tavg(h);
            l2(h)=lavg(h);
        end
        t2=t2./numTraj;
        l2=l2./numTraj;
        r=ct*kdt;
        s=cs*kds*sevprob*4/13;
        rds=r/s;
%         th_means(a)=exp(rds)*rds^(1-rds)*(gamma(rds)-gammainc(rds,rds));
        
        % th_dist=zeros(1,51);
        % th_dist_x=zeros(1,51);
        % for leng=1:length(th_dist)
        %      lengs=leng*cap-(cap-1);
        %      th_dist_x(leng)=lengs;
        %         th_dist(leng)=(lengs*s/r)*exp(-(s/r)*(lengs^2)/2);   
        % end 
            
    end
    cap=round(max(sst_lengths)/50);
    sst_lengths=sst_lengths*4/13/1000;
    p=histogram(sst_lengths(1:c-1),'Binwidth',.1,'Normalization','probability');
    area=0;
    
    
        th_dist_x=p.BinEdges+p.BinWidth/2;
        th_dist=zeros(1,length(th_dist_x));
        for i=1:length(th_dist_x)
             lengs=th_dist_x(i)*1000*13/4;
             th_dist(i)=(lengs*s/r)*exp(-(s/r)*(lengs^2)/2);   
        end 
    th_dist=th_dist./sum(th_dist);
    %need to scale theoretical curve because it is normalized as a
    %continuous distribution, but the histogram is normalized as discrete
    for i=1:p.NumBins
         coord=find(th_dist_x>=p.BinEdges(i),1);
         area=area+p.BinEdges(i)*th_dist(coord);
    end
    %th_dist=th_dist./sum(th_dist);
    hold on
    plot(th_dist_x,th_dist)
    title(titlen);
    xlabel('[T] (\muM)')
    ylabel('length (nm)')
    saveas(gca,'lengthdistns5.pdf');
end