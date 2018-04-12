%% PROGRAM INFORMATION
%---------------------------------------------------------------------
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date: Autumn 2017
%---------------------------------------------------------------------
%
% This file setup and starts the coupled 2.5-D model used for the spin-up
% test. It can be run with either 20 or 500 unsaturated columns. The
% program requires two setup files (for GRW and UZ) as well as some data
% specified in the "data" folder. More information about the coupling
% framwork can be found in runCoupled_V3p2.m
%
%---------------------------------------------------------------------
% LATEST EDITION
% Stripped version created for the publishing
%
%---------------------------------------------------------------------

% Number of x-pixels per zone in the GRW
np=50; % for the current problem: 50 --> 20 zones, 10 --> 500 zones


% Create the zonation map for the UZ columns
% Here: square boxes 
a=ones(np,np);
b=[];
for i=1:500/np
    b=[b;a];
    a=a*0+a(1)+1;
end
a=b;
b=[]; ix=max(max(a));
for j=1:100/np
    b=[b,a];
    a=a+ix;
end
LUVX=b; % and rename it

% Assinge landuse to two variables:
COUP.zoneMap=ones(500,100); % identifies differnt land uses (here: only 1)
COUP.zoneMapIndivid=LUVX; % identifies different zones which each has a UZ 
% column on top (smaller or equal to the land use map)

% Number of UZ cells
COUP.nrUZcells=max(max(LUVX));

% cell-index (y,x) for each of the COUP.nrUZcells
[x,y]=meshgrid(1:100,1:500);
COUP.UZindex=zeros(max(max(LUVX)),2);
for i=1:max(max(LUVX))
    ax=x(LUVX==i);
    ay=y(LUVX==i);
    %  LUV_v1(round(mean(ay)),round(mean(ax)))=8;
    COUP.UZindex(i,:)=[round(mean(ay)),round(mean(ax))];
end

clear CC ax ay c a i j x y
%% UZ and GRW model setups

% Get the basic setups from the setup files
% Groundwater flow model
[GRW.solver,GRW.grwms,GRW.boundary]=setupGRW_spinUp();
% Unsaturated zone model
[UZ.geom,UZ.param,UZ.solver,UZ.time,UZ.boundary,UZ.kin,UZ.alpha,UZ.n,UZ.Ssat,UZ.Sres,UZ.poro,UZ.aKR]...
    =setupUZ_spinUp(COUP.UZindex);

% Timing information
GRW.grwms.dt=3600*12;  % Time step (= coupling time step) [s]
GRW.grwms.dtStore=GRW.grwms.dt; % [s]
% Total run-time
GRW.grwms.tTot=3600*24*365*1; % [s] (here: 1 year)

% adjust so that we can run 50 years if needed
UZ.boundary.fluxTopMTRX=repmat(UZ.boundary.fluxTopMTRX,1,50);
UZ.boundary.changeTimes=0:3600:3600*24*365*50;

% Initial time for UZ models
UZ.time.dt=10; % [s]
UZ.time.dtMax=3600; % maximum time step size [s]

% How often should things be saved 
COUP.saveEach=7*2; % every X day
if np==10; COUP.saveEach=90*2; end


%% Initial conditions

% GRW from Steady State
load data/grwLevelSS_ISM grwlevelSS
grwlevelSS(grwlevelSS>49.9)=49.9; % adjust so that no flooding occurs
hGRW_Initial =grwlevelSS+GRW.grwms.z0; % create the initial head
clear hEnd


% Load the initial condition for the UZ
% Here: from steady state simulation of the individual columns and
% ranginging from bed-rock until land surface
load data/hUZ_Initial hUZ_Initial

% Are the unsaturated columns all homogeneous?
COUP.doUZhomo = false;

% Saving: pick name of directory to save in
sDir='SpinUp_run1';

reStartFrom=[];     % [] = start from the beginning

%% Run the 2.5-D model

runCoupled_V3p2(UZ,GRW,COUP,hUZ_Initial,hGRW_Initial,sDir,reStartFrom);
















