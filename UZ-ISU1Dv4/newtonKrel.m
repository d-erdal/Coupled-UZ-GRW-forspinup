% Part of a new version of the solution for the Richards equation
% Written May 2012 by Daniel Erdal
% Compacted version written May 2012

% Version 2.0
% Last Documented edition: 05.12.13 (DE)

% Calculation of the relative permeability.
% The resulting vector is size nz+2, meaning that both a fictive top and
% the bottom is in the vector.

function kr = newtonKrel(h,se,param,boundary)
% (input se is the effective saturation from the "calcSat" function)

% Initialize the resulting vector
kr=zeros(param.nzKR,1);


if strcmp(param.what,'VG')
    %% Calculate Krel with the Mualem Van Genucthen model
    al=param.alpha;
    n=param.n;
    aKR=param.aKR;
    
    a=al.*abs(h);
    kr(2:end-1)=((1-((a).^(n-1)).*se).^2).*(se.^aKR);
    kr(1)=boundary.kr; 
    
    
elseif  strcmp(param.what,'RG')
    %% Calculate Krel with the Russo-Gardner model 
    
    al=param.alpha;
    
    kr(2:end-1)=exp(-al.*abs(h));
    kr(1)=boundary.kr;    

   
end

% The bottom if free outflow:
% here we force an upsteam weighting
if boundary.doFreeDrain
    kr(1)=kr(2);
end
% The top: only used if the evaporation is broken off
% here we force an upsteam weighting that in this case (evaporation) is
% from below
kr(end)=kr(end-1);

% if positive pressure; full saturation; rel. perm. = 1
kr([boundary.BottomHead;h]>0)=1;
kr(kr>1)=1;


% limit the minimum value possible for the permeability
kr(kr<param.KrelMin)=param.KrelMin;



