%% PROGRAM INFORMATION
%---------------------------------------------------------------------
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date: August/September 2014,
%---------------------------------------------------------------------
%
% The groundwater solver solves the nonlinear equations for a 
% transient unconfined (phreatic) aquifer in horizontal 2D:
% n*(dh/dt) - div(K*(h-z)*grad(h)) = Q
% where n is porosity, K saturated conductivity, ground elevation, Q
% recharge and h groundwater level over NN.
%
% This is the setup file for the groundwater part of the 2.5-D model
%
%---------------------------------------------------------------------
% Output: solver - structure containing the solver setttings
%         grwms -   structure containing GRW-flow settings
%         boundary- structure containing GRW-boundary conditions%        
%
%---------------------------------------------------------------------
% LATEST EDITION
% Stripped version created for the publishing
%
%---------------------------------------------------------------------


function [solver,grwms,boundary]=setupGRW_spinUp()



%% Start file for running groundwater problem in 2D horizontal
% This start file is divded into X sections:


%% -----------Solver setup-------------------------------------

% How to solve nonlinear equations: 1) Newton 2) Fixpoint 3) Explicit GPU
solver.type =  1;
% For the fixpoint: use LU-factorization (see INFO in doConfined above)
solver.doLUfact = false;
% For the Newton method: dh for numerical jacobian
solver.delta_eps = 1e-6;
% Use a line search for the increment in the Newton method?
solver.doLineSearch = true;
% How to solve the linear equations: 1) MATLAB \ 2) Biconjugate gradient 3) gmres or GPU-Jacobi
% (not active if doLUfact and type==2 are true)
solver.how = 2;
% Print stuff?
solver.doPrint=true;
% Precondution with Explicit Euler?
solver.doExplicitPredcond = false;
% What maximun deltah is allwoed for explicit solver/precondutioners
solver.EulerMax=1e-5;

% Do steady state?
solver.doSS=false;

% Convergence for the absolute error
solver.cErr=1e-7;
% Maximum number of internal iterations
solver.maxJacobi=1000;
%---------------------------------------------


%% -----------Domain setup-------------------------------------

% INFO: The domain is built around a mask that contains ones and zeros for
% cells that are either in or outside the domain. Each cell on the boundary will
% be a boundary cell and will either have an assigned boundary (see below) or
% will bple defaulted to a no-flow. The mask can have any arbitrary shape that fits
% the grid. All later assigned fields will be cropped to the mask file and just have to
% have the right nx and ny.
boundary.mask=ones(500,100);

% Calculate ny and nx
nx=size(boundary.mask,2);
ny=size(boundary.mask,1);

grwms.ny=ny;
grwms.nx=nx;

% Discretization x and y (program runs best if dx=dy).
grwms.dx=10; % [L]
grwms.dy=10; % [L]

% Conductivity
load data/parameters_hetero Kms pinit poroBLOCK
Kms(pinit<0)=nan;
grwms.K=nanmean(Kms,3);

% Mean of two cells by: 1) Artihmetic, 2) Harmonic mean
grwms.Kmean=2;

% Porosity [-]
poroBLOCK(pinit<0)=nan;
grwms.poro=nanmean(poroBLOCK,3);

% Assing the bottom elevation [L]
load data/demV2 dem2

depth=50;
% lower the rivers from the DEM
dem2=dem2-min(min(dem2));
dem2wRivers=dem2;
dem2wRivers(:,1)=dem2(:,1)-1;
dem2wRivers(:,100)=dem2(:,100)-1;
dem2wRivers(1,:)=dem2(1,:)-1;

% Land-surface (only effects the rivers)
grwms.lsurf=dem2wRivers+depth;
% Bottom elevation
grwms.z0=grwms.lsurf-50;

%% -----------Time setup-------------------------------------
% THIS WILL BE OVERRIDDEN IN THE MAIN CODE LATER
% Total runtime [T]
grwms.tTot=[]; % seconds
% Basic time step [T]
grwms.dt=[];
% Time iterval for storing the output [T]
grwms.dtStore=[];
% Starting time
grwms.tStart=0;

% Do adaptive time stepping?
solver.doApdaptDT=false;

% DO STEADY STATE?-------
if solver.doSS
    grwms.poro=grwms.poro*0;
    grwms.dt=1;
    grwms.tTot=1;
    grwms.dtStore=1;
    solver.doExplicitPredcond = false;
    if solver.type>2
        error('No steady state for expicit solvers....')
    end
end
%-------------------------

%% -----------Boundary setup----------------------------------

% INFO: Boundary setup
% Two type of boundary conditions are implemented: Diriclet and Neumann
% The Recharge (Q) is Neumann, the other boundaries are either specified
% or defautls to a no-flow (Neumann with value 0).

% INFO: Temporal setup (all boundaries)
% Each boundary has a doX boolean that tells whether or not the variable should
% be updated in time (if not, only initial values are considered). The corresponding
% changeTX vector contains the time [T] after which a new boundary in the corresponding
% matrix/vector X is to be applied. Please note that the first value in changeTX HAS TO
% BE 0! If X is longer than changeTX, only the values corresponding to changeTX will be
% considered.

% INFO: Spatial setup (perimeter boundaries)
% The side boundaries are defiend using a perimeter matrix. Each cell that is a boundary
% cell can in the perimeter matrix be assigned a number. Each number is then associated
% with a certain boundary condition type and value. Each number can have a different value
% of the boundary, but temporal changes are assigned once for Diriclet and once for Neumann.
% Cells that are boundaries but have no assignment will later be assigned as no-flow.
% Please note: A cell that has more then one boundary can still only hold one boundary
% condition and this boundary condition is assigned to all boundaries! The only exeption is
% if the domain is a perfect box, then the edges can be diriclet and no flow.

%% Recharge boundary
% Here: recharge is a computed boundary from the UZ part of the coupled
% model; hence no assignment
 boundary.chngQ=false;
boundary.changeTQ=0;
boundary.QtempoScale=1;


%% Lateral boundaries
% Here: no latteral boundaries considered

% Assinge the perimeter cells an index
boundary.perimI=zeros(ny,nx);

% Which indecies are Diriclet boundaries
boundary.perimIdir=[];
% Which indecies are Neumann boundaries
boundary.perimIneu=[];

% Are there any internal diriclet boundaries
boundary.doInternal=false;

% Setup for Diciclet boundaries
boundary.chngDir=false;
boundary.changeTdir=0;
boundary.boundDirValues=[];

% Setup for Neumann boundaries
boundary.chngNeu=false;
boundary.changeTneu=0;%
boundary.boundNeuValues=[];


%% The river boundray
boundary.doRiver = true;

% create a riverMask (1 = precense of river)
boundary.riverMask=zeros(500,100);
boundary.riverMask(1:500,1)=1;
boundary.riverMask(1:500,end)=1;
boundary.riverMask(1,:)=1;

% Resoltion of the river mask 
% (can also be finer than the main grid!)
boundary.riverDX=grwms.dx;
boundary.riverDY=grwms.dy;

% Temporal setup of river boundary
boundary.chngRiver=false; % does it change over time?
boundary.changeTRiver=boundary.changeTdir;

% River parameters
boundary.Kriver=1e-6; % Conductivity of river bed
boundary.dRiver=0.05;    % Thickness of river bed

% Assign values for the river boundary 
% Either as fileds ny x nx or as scalars%
boundary.hRiver=0.05;

% RIVER ASSUMPTIONS
% The river head is measured from the land lsurf (hence, the hRiver should
% often be negative!)


%% -----------Wells setup-------------------------------------
%INFO: For each well, a location (y,x) need to be given together with
% start and stop times and pumpung rates (negative = extration). The program
% later assignes the well to a grid cell, so that the pumping rate hence is
% dependent on the dx-dy discretization.
% Please note that only one pumping rate is possible for each well, if more
% rates are required these will have to be assigend as new wells.

% Do we have wells
grwms.doWells=false;

% Number of wells
grwms.wells.nr=1;

% Position of wells [L]
grwms.wells.x = 500*2+2;
grwms.wells.y = 50*2+2;

% Start and Stop times for the wells [T]
grwms.wells.tStart = 0;
grwms.wells.tStop = inf;

% Pumping rates [L/T] ()
% 2015 edition: positive is extraction!!!!!!
grwms.wells.Q = 1e-4;


%% -----------Run the problem----------------------------------
% THIS IS ONLY A SETUP FILE; HENCE, WE DO NOT RUN THE PROBLEM





