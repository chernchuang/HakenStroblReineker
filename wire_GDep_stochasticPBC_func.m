function [MSD,D_num,pop,t]=wire_GDep_stochasticPBC_func(Geom,Para,Time)

if nargin<3
    Time.dt=0.1;
    Time.nt=1e2;
end
if nargin<2
    Para.nS=5e1;
    Para.GammaN=[0 0.1];
    Para.SigmaN=1;
    Para.w=10;
end
if nargin<1
    Geom.N=1e2;
end

% close all
% clear
% Time.dt=0.1;
% Time.nt=1e3;
% Para.nS=2e1;
% Para.GammaN=[0 0.1];
% Para.SigmaN=1;
% Para.w=10;
% Geom.N=2e2;

N=Geom.N;
nS=Para.nS;
GammaN=Para.GammaN;
SigmaN=Para.SigmaN;
NG=length(GammaN);

dt=Time.dt;
nt=Time.nt;

nlabel=-N/2+1:N/2;
t=(0:nt)*dt;

ind=toeplitz(1:N,[1 N:-1:2]);
iind=hankel(1:N,[N 1:N-1]);
iind=bsxfun(@plus,iind,(0:N-1)*N);

if isfield(Para,'w')&&(Para.w>0)
    phi0=exp(-nlabel.^2/Para.w^2)';phi0=phi0/norm(phi0);
else
    phi0=zeros(N,1);phi0((N)/2)=1;
end
if isfield(Para,'phi0')
    phi0=Para.phi0;
end

% Hs=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1); %OBC
Hs=diag(ones(N-1,1),1);Hs(1,N)=1; % PBC
Hs=Hs+Hs'; 
% expHs=expm(-1i*Hs*dt); % pure system propagator

phi0=phi0(ind); % PBC for the initial state

pop=zeros(N,nt+1,NG);
for kS=1:nS
    disp([kS nS])
    
    disorder=randn(N,1)*SigmaN;
    expHsd=expm(-1i*(Hs+diag(disorder))*dt); % pure system propagator
    
    for kG=1:NG
        phit=phi0;
        pop(:,1,kG)=pop(:,1,kG)+mean(reshape(abs(phit(iind)).^2,[N N]),2);
        for kt=1:nt
%             Gammat=sqrt(GammaN(kG))*randn(N,1);% stochastic term of HSR system-bath coupling
            Gammat=sqrt(GammaN(kG))*randn(N,N);% stochastic term of HSR system-bath coupling
            phit=bsxfun(@times,expHsd,exp(-1i*Gammat*dt))*phit;
%             phit=bsxfun(@times,expHs,exp(-1i*(Gammat+disorder)*dt))*phit;
            pop(:,kt+1,kG)=pop(:,kt+1,kG)+mean(reshape(abs(phit(iind)).^2,[N N]),2);
        end
    end
end

pop=pop/nS;

MSD=NaN(nt+1,NG);
for kG=1:NG
    MSD(:,kG)=nlabel.^2*pop(:,:,kG);
end
D_num=[zeros(1,NG);diff(MSD,1,1)/dt];

% 
% %%
% LW=2;
% Gcolor=get_lcolor(NG);
% 
% 
% figure
% yyaxis left
% c11=plot(t,MSD,'-','LineWidth',LW);
% [c11(:).Color]=deal(Gcolor{:});
% yyaxis right
% c12=plot(t,D_num,'--','LineWidth',LW);
% [c12(:).Color]=deal(Gcolor{:});
% yl12=get(gca,'ylim');
% set(gca,'ylim',[0 yl12(2)])
% 
% figure
% c2=semilogy(squeeze(pop(:,end,:)),'LineWidth',LW);
% [c2(:).Color]=deal(Gcolor{:});
% set(gca,'ylim',10.^[-10 1])
% tileF
