

close all
clear

w=0;
Time.dt=0.1;
Time.nt=1e3;
% Para.GammaN=[0 logspace(-2,2,9)];
Para.GammaN=0:1:10;
Para.nS=5e2;
Para.SigmaN=1;Geom.N=2e2;
Para.w=w;

[MSD,D_num,pop,t]=wire_GDep_stochasticPBC_func(Geom,Para,Time);Sto='_Stochastic';

% fn=sprintf('HSR_GDep_N%d_s%2.2e_w%2.2e%s.mat',Geom.N,Para.SigmaN,Para.w,Sto);
fn=sprintf('HSR_lGDep_N%d_s%2.2e_w%2.2e%s.mat',Geom.N,Para.SigmaN,Para.w,Sto);

if exist(fn,'file')
    a=load(fn);
    if numel(a.Para.nS)==1
        a.Para.nS=repmat(a.Para.nS,1,length(Para.GammaN));
    end
    nnS=Para.nS+a.Para.nS;
    MSD=(MSD.*Para.nS+a.MSD.*a.Para.nS)./nnS;
    D_num=(D_num.*Para.nS+a.D_num.*a.Para.nS)./nnS;
    pop=(pop.*permute(Para.nS,[1 3 2])+a.pop.*permute(a.Para.nS,[1 3 2]))./permute(nnS,[1 3 2]);
    Para.nS=nnS;
end

save(fn,'Time','Para','Geom','MSD','D_num','pop','t')
clear Time Para Geom MSD D_num pop t
