% Part of a new version of the solution for the Richards equation
% Written Dec 2011 by Daniel Erdal
% Compacted version written May 2012

% Version 2.0
% Last Documented edition: 21.12.16 (DE)--> ad-hoc 1 layer root uptake


% Evaluation of the Richards equation for the newton solver

% This means the evaluation of the following function:
% store = wc_new./param.poro.*param.storage.*(h-head_old);
% f = (store + wc_new-wc_old).*delzCS./dt...
%    -(K.plus.*((h(+)-h(x))./delzNN+1) - K.minus.*((h(x)-h(-))./delzNN+1));
%   -----
%   | + |
%   -----
%   | x |
%   -----
%   | 0 |
%   -----


function f=richardsNewton(h,kr,Ksat,storeX,delzNN,boundary)
nz=length(h);
dh=diff(h);

global doFlow

%% Weighting of permeabilities
% Ksat weighted harmonically before the solver starts and is an input
% (and is a vector nz+1 including the boundaries)
% Krel is weighted by upsteam weighting, meaning that the value on the
% upsteam node gets all the value.

if boundary.KrelWeight == 1
    % Weight by upsteam weighting of relative permeabilities
    % Create a vector of heads differences, with bottom and a top setup so no
    % incorrect weighting will occur on the boundaries (upstream is always
    % the cell that actually exists)
    dh2=ones(nz+1,1)*1e6;
    dh2(end)=boundary.ETbreak-h(end);  % as this is only used for ETbreak
    dh2(2:end-1)=dh;
    
    if ~boundary.doFreeDrain
        dh2(1)=h(1)-boundary.BottomHead;
    end
    
    % if hpu(x) is true; use the upper node: Krel(x)=kr(x)
    % if hpu(x) is false; use the lower node: Krel(x)=kr(x+1)
    hpu=dh2>-delzNN;
    
    b=2:nz+2;
    Krel=kr(1:end-1);
    Krel(hpu)=kr(b(hpu));
elseif boundary.KrelWeight == 2
    % Do arithmetic averaging
    Krel=0.5*(kr(1:end-1)+kr(2:end));
    
elseif boundary.KrelWeight == 3
    % Do harminic averaging
    Krel=(2*kr(1:end-1).*kr(2:end))./(kr(1:end-1)+kr(2:end));
    
elseif boundary.KrelWeight == 4
    % geometric averaging
    Krel=sqrt(kr(1:end-1).*kr(2:end));
    
end


% Compute the output
% (one output, since upper for one is lower for the other)
K=Ksat.*Krel;


%% Internal nodes

% initialze the result vector
f=zeros(nz,1);

% calculate the fluxes
khx=K(2:end-1).*(dh./delzNN(2:end-1)+1);
% calculate the divergence of the fluxes
dflux=diff(khx);
% calculate the Richards equation
f(2:end-1)=storeX(2:end-1)-dflux;


%% The boundaries
% Top: flux
 if boundary.TopFlux>0 % if evaporation
    % chose the flux boundary or the ETbreak pressure head:    
    fbnd=max((K(end).*((boundary.ETbreak-h(end))./delzNN(end)+1)),-boundary.TopFlux);
    
    f(end)=storeX(end)-(fbnd-khx(end));
 else
      f(end)=storeX(end)-(-boundary.TopFlux-khx(end));
 end

 if boundary.doAdHoc1LayerRoots
     % a very odd one layer root water uptake. Done to match how fluxes can
     % be limited in parflow.
     % REMOVE WHEN NO LONGER IN USE
     if h(end)>-100
         rwu=boundary.adHocUptake;
     elseif h(end) >-150
         % feddes limiter
         rwu=(0.95*(h(end)+150)/50+0.05)*boundary.adHocUptake;
     else
         % still keep somethings, but keep it small
         rwu=0.05*boundary.adHocUptake;
     end   
     
    f(end)=f(end)+rwu; 
     
 end
 

% Bottom: head or free drain
if boundary.doFreeDrain   
    
    if boundary.freeDrainOption == 1
         % dh/dz=0 at the bottom, hence only the K term left
         % Also known as the zero order extrapolation or simply as gravity
         % driven outflow
        f(1) = storeX(1)-(khx(1)-K(1));        
        
    elseif boundary.freeDrainOption == 2
        % Alt 2, dh/dz is the same at both faces of the last cell
        % Also known as a 1st order extrapolation of the heads in the last
        % 2 cells onto a ghost cell
        f(1) = storeX(1)-(khx(1)-K(1)*(khx(1)/K(2)));%           
        
    elseif boundary.freeDrainOption == 3
         % The bottom is a gravel layer and the outflow is either gravity
         % (option 1 above) if the bottom layer is fully saturated or no 
         % outflow at all if saturation < 1
         
         if strcmp(boundary.pWhat,'VG')
             % Van Genuchten
             s1=1./((1+(boundary.alphaB*abs(h(1))).^boundary.nB).^boundary.mB);
         elseif strcmp(boundary.pWhat,'RG')
             % Russo-Gardner          
             s1=(exp(-0.5*boundary.alphaB*abs(h(1))).*(1+0.5*boundary.alphaB*abs(h(1)))).^(2/(2+boundary.aKR));
         end
         
         s1(s1>1)=1;
         s1(h(1)>0)=1;
         
         if s1>=0.9 && doFlow
             % if saturation above 0.9: full drain%         
             x=interp1([0.9 0.93 0.95 0.97 1],[0 0.1 0.5 0.9 1],s1,'PCHIP');
             %x=interp1([0.9 0.93 0.95 0.97 1],[0 0.1 0.5 0.9 1],s1,'cubic');
            f(1) = storeX(1)-(khx(1)-K(1)*x); 
         else
             % no flow bottom
             f(1) = storeX(1)-(khx(1)-0);  
         end       
        
    end
else
       % if no free drain, then it is a head bottom   
    if boundary.doBottomFlux0
        f(1) = storeX(1)-(khx(1)-0);
        
    elseif boundary.noBackFlow        
        % for the lysimeter cases, make sure no flow goes back into the system
        fbnd=max(0,K(1).*((h(1)-boundary.BottomHead)./delzNN(1)+1));
        f(1) = storeX(1)-(khx(1)-fbnd);        
        
    else
        % use the specified head boundary
      f(1) = storeX(1)-(khx(1)-K(1).*((h(1)-boundary.BottomHead)./delzNN(1)+1));
%     
    % temporary test: ghost cell approach to bottom boundary
%     f(1) = storeX(1)-(khx(1)-Ksat(1).*((h(1)-(boundary.BottomHead+delzNN(1)))./(2*delzNN(1))+1));%+boundary.adHocBottomFlux;
    end
end


f=f+boundary.adHocBottomFlux;

 %  f(end)=storeX(end)-(-boundary.TopFlux-khx(end));

