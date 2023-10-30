% plots main figure for figure 3d, which is the axes

%hold on
MaxT=10000000;
numTraj=200;
%severing repair parameters
ks=.54;
kt=.11;
x0=1;
L=25;
cts=.01:1:5.01; %subunits/s  %tubulin concentrations
cat=.6;  %.004;%/s %catastrophe rate
res=.7;  %.02; %rescue rate
vs= 1300;%152; %shrinkage rate in catastrophe
cs=.05; %severing protein concentration
kds=12/6000; %severing rate
kdt=340; % growth rate

%Compute severing probability and mean severing time
pt=kt*cts./(kt*cts+ks);
ps=(ks./(kt*cts+ks));
tstep=zeros(size(cts));
for a=1:length(cts)
    tstep(a)=1/(ks+kt*cts(a));
end
r=pt./ps;
one_vec=1;
sevprobs=((pt./ps).^(x0)-1)./( ((pt./ps).^L)-1); %probability
means1=tstep.*(((r+1)./(r-1)).*(((r.^L+1)./(r.^L-1)).*L-((r.^x0+1)./(r.^x0-1)).*x0));

%initialize variables
sevdata=0;
means=zeros(1,length(cts));
stdevs=zeros(1,length(cts));
sst=800; %time at which steady state is reached
stop_time=12000;
for n=[2,3]
        if n==2
            comp=true;
            titlen='with repair,w/o delay';
        elseif n==3
            comp=false;
            titlen='w/o repair';
        end
    for a=1:length(cts)
        mean1=means1(a);
        prob=zeros(1,3000000);
        sevprob=sevprobs(a);
        ct=cts(a);
        lengths=zeros(1,numTraj);
        c=1;
        sst_lengths=zeros(1,(stop_time-sst+1)*numTraj*20);
        cat_cut=0;
        rescs=0;
        for j=1:numTraj
            l=zeros(1, MaxT);
            t=zeros(1, MaxT);
            dt=zeros(1, MaxT);
            severing=zeros(1,length(prob));
            numsev=0;
            resc=true;
            cata=false;
                sev=0;

            for i=1:MaxT
                seva=0; %reset amount severed to 0 for each step
                 k2=cs*(l(i)-1)*4/13*kds; %multiply by length of mt to multiply likelihood of landing by that of severing
                k1=ct*kdt; %growth
                k4=vs;%shrink
                k3=cat;
                k5=res;
                if l(i)<=1
                    k2=0;
                end
                if l(i)==0
                    k4=0;
                end
                if resc
                    k0=k1+k2+k3;
                else
                    k0=k2+k4+k5;
                end
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
                if p<= (k2/k0)
                    l(i+1)=l(i); %length remains until a change occurs- in this time step only a severing protein lands
                    sevrand=rand;
                    if comp==false
                        sevprob=1;
                    end
                    if sevrand<sevprob
                        w=(k2/k0)/(l(i)-1);  %break up remaining possible p into segments
                        c1=0; %lowest possible random variable
                        c2=c1+w; %increment for each possible segment
                        for len=1:l(i) %for every possible length
                            s=len;
                            if p>c1 && p<=c2 %if p is in corresponding interval, land at this position
                                if comp==false
                                    l(i+1)=s; 
                                    sevdata=sevdata+1;
                                    numsev=numsev+1;
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
                            severing(sev)=t(i+1);
                        end
                    end
                elseif resc && p<=(k1+k2)/k0
                    l(i+1)=l(i)+1;
                elseif resc
                    resc=false;
                    cata=true;
                    l(i+1)=l(i);
                elseif cata && p<=(k4+k2)/k0
                    l(i+1)=l(i)-1;
                elseif cata %chance is k5
                    resc=true;
                    cata=false;
                    if(l(i)>0)
                    cat_cut=cat_cut+l(i);
                    rescs=rescs+1;

                    end
                    l(i+1)=l(i);
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
                            sevdata=sevdata+1;
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
        %plot lvt trajectories
        if n==2
            figure(1)
            plot(t(1:i),l(1:i),'DisplayName',num2str(ct))
            hold on
        end
        if n==3
            figure(4)
            plot(t(1:i),l(1:i),'DisplayName',num2str(ct))
            hold on
        end
    end
    %figure

    %plots the means
    figure(2)
    %load("dyninstinvivomean.fig")
    %hold on
    plot(cts,means,'-o','DisplayName',titlen)
    %histogram(sst_lengths(1:c-1))
    hold on

    %plot CV
    figure(3)
    %load("dyninstinvivocv.fig")
    %hold on
    plot(cts,stdevs./means,'-o','DisplayName',titlen)
    hold on

end

% add axis labels and legends
figure(1)
title('with repair')
legend;

figure(2)
xlabel('[T] (\muM)');
ylabel('mean steady state length (nm)');
legend;

figure(3)
xlabel('[T] (\muM)');
ylabel('CV (stdev/mean)');
legend;

figure (4)
title('without repair)')
legend;


save('dyn_inst_ct.01_sstlengths','sst_lengths')