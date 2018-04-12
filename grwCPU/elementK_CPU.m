% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% For further referecens see the main file runGRWunconf.m



function [kii,kit,kib,kil,kir]=elementK_CPU(ki,kt,kb,kl,kr,dx,dy)

% Harmonic avarage of the conductivies
kit=(2.*kt.*ki./(kt+ki))./(dy.^2);      kit(isnan(kit))=0;
kib=(2.*kb.*ki./(kb+ki))./(dy.^2);      kib(isnan(kib))=0; 
kil=(2.*kl.*ki./(kl+ki))./(dx.^2);      kil(isnan(kil))=0; 
kir=(2.*kr.*ki./(kr+ki))./(dx.^2);      kir(isnan(kir))=0; 

% for the diriclet boundaries
kii=ki./(dx./2.*dy);

