% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% For further referecens see the main file runGRWunconf.m

function a=forNorm_CPU(ji,jt,jb,jl,jr,dxi,dxt,dxb,dxl,dxr,f0)

aa=f0+ji.*dxi+jt.*dxt+jb.*dxb+jl.*dxl+jr.*dxr;
a=abs(aa).^2;









