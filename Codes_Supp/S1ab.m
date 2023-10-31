%least squares but with theory
cts=[0,2,5];
Ls=24:30; %Nd
kts=.05:.01:.2; %k_T in paper
kss=.4:.01:.7; %k_r in paper
x0s=1; %start
lstsqrsm=zeros(length(kss),length(kts),length(Ls),length(x0s));
for w=1:length(x0s)
    x0=x0s(w);
    for z=1:length(Ls)
        L=Ls(z);
        for x=1:length(kts)
            kt=kts(x);
            for y=1:length(kss)
                ks=kss(y);
                kT=kt.*cts;
                pt= kT./(kT+ks);
                ps= ks./(kT+ks);
                r=pt./ps;
                one_vec=ones(size(cts));
                k=ones(size(cts)).*x0;
                n=ones(size(cts)).*L;
                tstep=zeros(size(cts));
                %syms j
                for a=1:length(kT)
                    tstep(a)=1/(ks+kT(a));
                end

                C5=tstep.*(((r+one_vec)./(r-one_vec)).*(((r.^n+one_vec)./(r.^n-one_vec)).*n-((r.^k+one_vec)./(r.^k-one_vec)).*k));
                for i=1:length(kT)
                    if kT(i)==ks
                        C5(i)=(L^2-x0^2)/3;
                    end
                end
                data=[32.3,9.2; %0: mean, std 
                      75.4,30.8; %2: mean, std 
                      193.8,101.5]; %5: mean, std 
                lstsqrm=sum(((data(:,1))'-C5).^2);
                lstsqrsm(y,x,z,w)=lstsqrm;
            end
        end
    end
end
% plot heatmaps of values
kt_ind=7; %index for value of kt you want to plot wiht
figure
lstsqrsm(isnan(lstsqrsm))=1000;
lstsqrsm(lstsqrsm>1000)=1000;
imagesc('XData',Ls,'YData',kss,'CData',squeeze(lstsqrsm(:,kt_ind,:,1)));

colorbar
xlabel('N_D')
ylabel('k_r')
title("k_t="+num2str(kts(kt_ind)))
%%
Lind=2; %index for value of Nd you want to plot with
figure
imagesc('XData',kts,'YData',kss,'CData',squeeze(lstsqrsm(:,:,Lind,1)));
colorbar
xlabel('k_t')
ylabel('k_r')
title("N_D="+num2str(Ls(Lind)))

