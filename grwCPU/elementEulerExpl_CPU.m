% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% For further referecens see the main file runGRWunconf.m


function [flowB,dtMinX]=elementEulerExpl_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,zt,zb,zl,zr,Dbound,ki,zi,inflow,poro,kRi,Emin)
mwt=0.01; % minimum water table
%cflMin=0.8;
%dx=10; dy=10;
qt=kt.*((hi+ht)./2-zt).*(ht-hi);
ix=(hi+ht)./2-zt<mwt;
if any(ix)
    qt(ix)=kt(ix).*mwt.*(ht(ix)-hi(ix));
end

qb=kb.*((hi+hb)./2-zb).*(hb-hi); 
ix = (hi+hb)./2-zb<mwt;
if any(ix)
    qb(ix)=kb(ix).*mwt.*(hb(ix)-hi(ix));   
end


ql=kl.*((hi+hl)./2-zl).*(hl-hi); 
ix = (hi+hl)./2-zl<mwt;
if any(ix)
    ql(ix)=kl(ix).*mwt.*(hl(ix)-hi(ix));   
end

qr=kr.*((hi+hr)./2-zr).*(hr-hi);
ix = (hi+hr)./2-zr<mwt;
if any(ix)
    qr(ix)=kr(ix).*mwt.*(hr(ix)-hi(ix));  
end



b0=zeros(size(hi));
ix = ~isnan(Dbound);
if any(ix)
    % the boundaries are Diriclet
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....
    b0=ki.*(Dbound-zi).*(Dbound-hi);     
    ix2 = Dbound-zi<mwt;
    if any(ix2)
        b0(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-hi(ix2));  
    end
end


flowB=(qt+qb+ql+qr+b0+inflow-kRi.*hi)./poro;

%%

% find min and max of h
% hmax=max(hi,max(ht,max(hb,max(hl,hr))));
% hmin=min(hi,min(ht,min(hb,min(hl,hr))));

dtMinX=Emin./abs(flowB);
























