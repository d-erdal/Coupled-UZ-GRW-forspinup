%% PROGRAM INFORMATION
%---------------------------------------------------------------------
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date: 2011-2012,
%---------------------------------------------------------------------
%
% The unsaturated solver solves the mixed form of the Richards equation for
% a transient 1-D column subject to a range of different transient boundary
% conditions.
%
% This is the setup file for the unsaturated part of the 2.5-D model
%
%---------------------------------------------------------------------
% Input: cellindx - cell with the x-y location of each of the UZ-models 
%
% Output: geom   - structure containing UZ-geometry settings
%         param  - structure containing UZ-parameters and parameter info
%         solver - structure containing the solver setttings
%         time   - structure containing the time info
%         boundary - structure containing UZ-boundary conditions
%         kin    - matrix nrUZ-cells x nrUZ-models with conductivity values
%         alpha  - same as above but for alpha
%         n      - same as above but for n
%         Ssat   - same as above but for Ssat
%         Sres   - same as above but for Sres
%         poro   - same as above but for porostiy
%         aKR    - same as above but for turtuosity factor
%
%---------------------------------------------------------------------
% LATEST EDITION
% Stripped version created for the publishing
%
%---------------------------------------------------------------------



function [geom,param,solver,time,boundary,kin,alpha,n,Ssat,Sres,poro,aKR]=setupUZ_spinUp(cellIndx)

%% -------------- Solver options------------------------------------------- 
% Do the line search?
solver.doLineSearch = true;

% Convergence criterium for the iterations
solver.eps = 1.0E-7;

% Delta h for the finite difference of the retention function [m]
% (NR and P-N)
solver.epsh = 1.0E-7;

% Should the solver be used with EnKF or MCMC?
solver.doEnKF=false;
solver.doDream=false;

% Print run-time information
solver.doPrint=true;

% Convergence check 1: reduce the time step
solver.convBreak1=300;
% Convergence check 2: break out
solver.convBreak2=1500;

% if the jacobian is ill conditioned/ singular should the program break out
% (true) or reduce the time step and try again (false)
solver.jacSingBreak=false;

% ------------------------------------------------------------------------ 


%% -------------- Geometry parameters -------------------------------------

% Number of layers in problem (spatial resolution)
geom.nz= 50;

% Height of the domain H [m]
geom.heighcell=50;

% Grid: delzCS is the size of each of the cells in the grid
% (sum(geom.delzCS)=geom.heighcell and length(geom.delzCS)=geom.nz)
load data/dzScaleCLMtop dzX
geom.delzCS=dzX;

if size(geom.delzCS,2)~=1
    geom.delzCS=geom.delzCS';
end

% Grid: node to node distance
% This is automated based on the delzCS and need not be altered
% (sum(geom.delzNN)=geom.heighcell and length(geom.delzCS)=geom.nz+1)
geom.delzNN=zeros(geom.nz+1,1);
geom.delzNN(1)=geom.delzCS(1)/2;
for i=2:length(geom.delzCS)
    geom.delzNN(i)=(geom.delzCS(i-1)+geom.delzCS(i))/2;
end
geom.delzNN(end)=geom.delzCS(end)/2;

% ------------------------------------------------------------------------ 


%% -------------- Parmeterization parameters ------------------------------

param.what='VG';
% van Genuchten model; VG
    %(param. alpha, n, Ssat, Sres, aKR, poro)
% Russo-Gardner: RG
    %(param. alpha, Ssat, Sres, aKR, poro)

load data/parameters_hetero Kms poroBLOCK 
    
% Hydraulic conductivity field kf [m/s].
kin=zeros(geom.nz,size(cellIndx,1));
for i=1:size(cellIndx,1)
    kin(:,i) = squeeze(Kms(cellIndx(i,1),cellIndx(i,2),:));
end

% Porosity [-]
poro=zeros(geom.nz,size(cellIndx,1));
for i=1:size(cellIndx,1)
    poro(:,i) = squeeze(poroBLOCK(cellIndx(i,1),cellIndx(i,2),:));
end

% Inverse entry pressure parameter alpha [1/m]
alpha=zeros(geom.nz,size(cellIndx,1));
alpha(poro==0.1)=2.0; % rock
alpha(poro==0.465)=1.28; % soil
alpha(poro==0.43)=14.0; % gravel


% Shape factor for van Genuchten model
n=zeros(geom.nz,size(cellIndx,1));
n(poro==0.1)=2.1; % rock
n(poro==0.465)=1.6; % soil
n(poro==0.43)=2.6; % gravel
% param.m=1-1./param.n;

% Saturated saturation and residual saturation
Ssat=1.0 * ones(geom.nz,size(cellIndx,1));
Sres=0.01 * ones(geom.nz,size(cellIndx,1));

% Turtousity parameter (also known as tau)
aKR= 0.5 * ones(geom.nz,size(cellIndx,1));

% Storage parameter for saturated conditions
param.storage = 1.0E-6;

% How are the relative permeabilities weighted?
% 1: upsteam
% 2: arithmetic
% 3: harmonic
% 4: geometric
param.KrelWeight=1;

% What minimum value can the relative permeability take?
% (param.KrelMin=0 --> no effect)
param.KrelMin=0;

% ------------------------------------------------------------------------ 


%% -------------- Boundry conditions --------------------------------------

% head boundary condition h at the bottom [m]
boundary.headBottom=0; 
boundary.doBottomFlux = false;
% Use a free drainage bottom boundary? (means the boundary.headBottom is
% meaningless)
boundary.doFreeDrain=false;
% what free boundary? 
% 1: gravity drain / 0-order extrapolation
% 2: equal flux / 1st-order extrapolation
% 3: on/off boundary condition
boundary.freeDrainOption = 1;

% flux boundary at the top [m/s]. Negative sign: Downwards direction
load data/PandET2009mps Ptot ETtot
% Scale precipitiation (Ptot) and ET (ETtot) to give the net flux
boundary.fluxTopMTRX(1,:)=Ptot*0.82+ETtot*0.53; % --> 254 mm/year 
boundary.adHocTopMTRX(1,:)=ETtot*0.0; 

% and time when to change it. First item in boundry.changeTimes
% must always be 0, and if only 0, then static boundry condition
%cumsum(ones(1,length(boundary.fluxTop)-1)*3600);   % change every hour
boundary.changeTimes = [0,cumsum(ones(1,length(Ptot)-1)*3600)];

% when should the evopration flux be replaced with a head boundary 
%(boundary.ETbreak = -inf--> never change)
boundary.ETbreak=-20000; % [m]

% if using a head boundary, should water be stopped from flowing back into
% the system? (typically if simulation a pressure bottomed lysimeter)
boundary.noBackFlow=false;

% Initial head condition: pressure head
% boundary.headinit=-0.5*ones(geom.nz,1);
boundary.doBottomFlux0=false;

% Set a few non-used features to 0 / false
boundary.adHocBottomFlux=0;
boundary.doAdHoc1LayerRoots=false;

% ------------------------------------------------------------------------ 

%% -------------- Time setup ---------------------------------------------

% % starting time of simulation [s]
% time.tStart=0;
% 
% % total time of simulation [s]
% time.tStop = 365*24*3600;  % run for one year

% initial time step; on or one for each boundary [s]
time.dt =3600;

% Result output 
% tStore can be a number (store every tStor seconds) 
% or a vector of exact saving times
time.tStore = 3600;

% Use automatic time adaptation?
time.doAdaptiveTime = true;
% maxinum time step size [s]
time.dtMax=3600;
% minimum time step size [s]
time.dtMin=1;
 
% Maximum (over the full simulation time) number of iterations allowed
% (time.MaxIter=inf --> no maximum iteration limit)
time.MaxIter=inf;
% ------------------------------------------------------------------------

%% ----------------- Start the solver ------------------------------------
 
% THIS IN ONLY A SETUP FILE, SO JUST EXIT HERE





