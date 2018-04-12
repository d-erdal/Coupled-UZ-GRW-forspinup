% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% For further referecens see the main file runGRWunconf.m


function [f0]=elementf0_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,zt,zb,zl,zr,poro,Dbound,ki,zi,hOld,inflow,dt,kRi)
mwt=0.01; % minimum water table

at0=kt.*((hi+ht)./2-zt).*(ht-hi); 
ix=(hi+ht)./2-zt<mwt;
at0(ix)=kt(ix).*mwt.*(ht(ix)-hi(ix));     


ab0=kb.*((hi+hb)./2-zb).*(hb-hi); 
ix=(hi+hb)./2-zb<mwt;
ab0(ix)=kb(ix).*mwt.*(hb(ix)-hi(ix));  

al0=kl.*((hi+hl)./2-zl).*(hl-hi);   
ix=(hi+hl)./2-zl<mwt;
al0(ix)=kl(ix).*mwt.*(hl(ix)-hi(ix));   

ar0=kr.*((hi+hr)./2-zr).*(hr-hi); 
ix=(hi+hr)./2-zr<mwt;
ar0(ix)=kr(ix).*mwt.*(hr(ix)-hi(ix));  

b0=zeros(size(ki));
ix=~isnan(Dbound);
if any(ix)
    % the boundaries are Diriclet
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....
    b0(ix)=ki(ix).*(Dbound(ix)-zi(ix)).*(Dbound(ix)-hi(ix));
    
    ix2=Dbound-zi<mwt;
    b0(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-hi(ix2));       
   
end


% functions
f0=poro.*(hi-hOld)-dt.*(at0+ab0+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;





