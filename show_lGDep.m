
clear
close all

N=200;SigmaN=1;Sto='_Stochastic';
nlabel=-N/2+1:N/2;

fn=sprintf('HSR_lGDep_N%d_s%2.2e_w%2.2e%s.mat',N,SigmaN,0,Sto);
load(fn)

%% delta function, Gamma dependence
LW=3;MS=12;FS=40;

Glist=Para.GammaN;
NG=length(Glist);
DG=Glist(2)-Glist(1);
Gcolor=get_lcolor(NG,1.5,0.8);

km=20; % moving average parameter

kts=[50 180];
% kts=[50 200];

Deltat=3*Time.dt; % range of time stamps being highlighted
yticks1=0:3:10;
yl=[0 2.2];
f1=figure;
c1=plot([0 movmean(t,km)],[zeros(1,NG);movmean(D_num,km,1)]/2,'-','LineWidth',LW);
[c1(:).Color]=deal(Gcolor{:});
a1=gca;
set(a1,'xlim',[0 20],'ylim',yl)
set(a1,'FontSize',FS,'TickLabelInterpreter','LaTeX','LineWidth',2)
xlabel('Time $Jt$','FontSize',FS,'Interpreter','LaTeX')
ylabel('$D(t)/Ja^2$','FontSize',FS,'Interpreter','LaTeX')
a2=axes('Position',[0.7 0.65 0.05 0.2]);
for kG=1:NG
    patch('Vertices',[0 0;0 1;1 1;1 0].*[1 DG]+[0 kG-1]*DG,'Faces',1:4,'FaceColor',Gcolor{kG},'EdgeColor','none')
end
set(a2,'xtick',[],'ytick',yticks1,'FontSize',FS/4*3,'TickLabelInterpreter','LaTeX','LineWidth',2)
ylabel('$\Gamma/J$','FontSize',FS/4*3,'Interpreter','LaTeX')
box on

xticks2=0:2:10;
f2=figure;
hold on
c2=plot(Para.GammaN,D_num(kts,:)/2,'o-','LineWidth',LW,'MarkerSize',MS);
hold off
box on
set(gca,'xlim',[-1 11])%,'xtick',10.^(-2:2),'ylim',yl)
set(gca,'FontSize',FS,'TickLabelInterpreter','LaTeX','LineWidth',2)
set(gca,'xtick',xticks2)
xlabel('$\Gamma/J$','FontSize',FS,'Interpreter','LaTeX')
ylabel('$D(t'')/Ja^2$','FontSize',FS,'Interpreter','LaTeX')
h2=legend('Short Time','Long Time');
set(h2,'FontSize',FS,'Interpreter','LaTeX','LineWidth',2)

figure(f1);
for kt=1:length(kts)
    patch(a1,'Vertices',[t(kts(kt))-Deltat yl(1);...
                         t(kts(kt))+Deltat yl(1);...
                         t(kts(kt))+Deltat yl(2);...
                         t(kts(kt))-Deltat yl(2)],'Faces',1:4,...
             'EdgeColor','none','FaceColor',c2(kt).Color,'FaceAlpha',0.5)
end
figure(f2);

tileF([1 2])



