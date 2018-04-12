%% PROGRAM INFORMATION
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% For further referecens see the main file runGRWunconf.m
% Creations date:  March./April 2015


function [ji,jt,jb,jl,jr,f0]=elementJacobian_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,zt,zb,zl,zr,poro,Dbound,ki,zi,hOld,inflow,dt,kRi)
mwt=0.01; % minimum water table
e=1e-10;

% Inflow from the top cell==============================
at0=kt.*((hi+ht)./2-zt).*(ht-hi);
atT=kt.*((hi+ht+e)./2-zt).*(ht+e-hi);
atI=kt.*((hi+e+ht)./2-zt).*(ht-(hi+e));
ix=(hi+ht)./2-zt<mwt ;
if any(ix)
    % If the water table too low that (h-z) is fixed
    at0(ix)=kt(ix).*mwt.*(ht(ix)-hi(ix));
    atT(ix)=kt(ix).*mwt.*(ht(ix)+e-hi(ix));
    atI(ix)=kt(ix).*mwt.*(ht(ix)-(hi(ix)+e));
end

% Inflow from the bottom cell==============================
ab0=kb.*((hi+hb)./2-zb).*(hb-hi);
abB=kb.*((hi+hb+e)./2-zb).*(hb+e-hi);
abI=kb.*((hi+e+hb)./2-zb).*(hb-(hi+e));
ix=(hi+hb)./2-zb<mwt;

if any(ix)
    ab0(ix)=kb(ix).*mwt.*(hb(ix)-hi(ix));
    abB(ix)=kb(ix).*mwt.*(hb(ix)+e-hi(ix));
    abI(ix)=kb(ix).*mwt.*(hb(ix)-(hi(ix)+e));
end

% Inflow from the left cell==============================
al0=kl.*((hi+hl)./2-zl).*(hl-hi);
alL=kl.*((hi+hl+e)./2-zl).*(hl+e-hi);
alI=kl.*((hi+e+hl)./2-zl).*(hl-(hi+e));

ix=(hi+hl)./2-zl<mwt;

if any(ix)
    al0(ix)=kl(ix).*mwt.*(hl(ix)-hi(ix));
    alL(ix)=kl(ix).*mwt.*(hl(ix)+e-hi(ix));
    alI(ix)=kl(ix).*mwt.*(hl(ix)-(hi(ix)+e));
end

% Inflow from the right cell==============================
ar0=kr.*((hi+hr)./2-zr).*(hr-hi);
arR=kr.*((hi+hr+e)./2-zr).*(hr+e-hi);
arI=kr.*((hi+e+hr)./2-zr).*(hr-(hi+e));

ix = (hi+hr)./2-zr<mwt;
if any(ix)
    ar0(ix)=kr(ix).*mwt.*(hr(ix)-hi(ix));
    arR(ix)=kr(ix).*mwt.*(hr(ix)+e-hi(ix));
    arI(ix)=kr(ix).*mwt.*(hr(ix)-(hi(ix)+e));
    
end

% Inflow from the boundaries==============================
b0=zeros(size(ki)); bI=b0;
ix = ~isnan(Dbound);
if any(ix)
    % the boundaries are Diriclet
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....
     b0(ix)=ki(ix).*(Dbound(ix)-zi(ix)).*(Dbound(ix)-hi(ix));
     bI(ix)=ki(ix).*(Dbound(ix)-zi(ix)).*(Dbound(ix)-(hi(ix)+e));
       
    ix2=Dbound-zi<mwt;
    
    if any(ix2)
        b0(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-hi(ix2));
        bI(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-(hi(ix2)+e));
    end
end


% Newton-Rahpson functions
f0=poro.*(hi-hOld)-dt.*(at0+ab0+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fI=poro.*(hi+e-hOld)-dt.*(atI+abI+alI+arI+bI)-dt.*inflow+dt.*kRi.*(hi+e);

fT=poro.*(hi-hOld)-dt.*(atT+ab0+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fB=poro.*(hi-hOld)-dt.*(at0+abB+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fL=poro.*(hi-hOld)-dt.*(at0+ab0+alL+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fR=poro.*(hi-hOld)-dt.*(at0+ab0+al0+arR+b0)-dt.*inflow+dt.*kRi.*hi;

% Construct the Jacobian entries 
% main diagonal
ji=(fI-f0)./e;
% first positive off diagonal
jt=(fT-f0)./e;
% first negative off diagonal
jb=(fB-f0)./e;
% second negative off diagonal
jl=(fL-f0)./e;
% second positive off diagon
jr=(fR-f0)./e;





