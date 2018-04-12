% Part of a new version of the solution for the Richards equation
% Written Dec 2011 by Daniel Erdal
% Compacted version written May 2012

% Version 2.0
% Last Documented edition: 30.05.12 (DE)

% Ouput is saturation (saft, for calculating water content)
% and effective saturation (se, for the rel. perm. function)

function [satf,se]=calcSat(h,param)


if strcmp(param.what,'VG')
    %% Calculate saturation according to the Van Genuchten model
      
    
    % first calculation
    a=param.alpha.*abs(h);
    
    % effective saturation
    satf=1./((1+a.^param.n).^param.m);
    satf(h>0)=1;
    satf(satf>1)=1;
    se=satf;
    
    % saturation
     satf=satf.*param.sDiff+param.Sres;
    
elseif  strcmp(param.what,'RG')
    %% calculate saturation accoring to the Russo-Gardner function
    
       
    a=0.5.*param.alpha.*abs(h);
    
    satf=(exp(-a).*(1+a)).^(2./(2+param.aKR));
    
    satf(h>0)=1;
    satf(satf>1)=1;
    se=satf;
    
    % saturation
    satf=satf.*param.sDiff+param.Sres;
    
    
    
end


