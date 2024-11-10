function lcc=get_lcolor(N,k,fac)
if nargin<2
    k=10;
end
if nargin<3
    fac=0.8;
end
lcolor=(jet(N)*k+1)/(k+1);
lcolor=lcolor*fac;
lcc=mat2cell(lcolor,ones(N,1),3);

