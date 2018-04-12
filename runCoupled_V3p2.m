%% PROGRAM INFORMATION
%---------------------------------------------------------------------
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date: Spring 2017
%---------------------------------------------------------------------
%
% The program combines a 2-D phreatic aquifer model (in the following: GRW) 
% with a small collection of 1-D unsaturated zone (Richards equation) models 
% (in the following: UZ). The UZ-models ranges from the groundwater table
% up to the land surface. First, the UZ models are run to aquire a recharge 
% that is then reginonlized to the full lateral model domain an used as a top
% boundary for the GRW-model. As a second step, the UZ model domains are
% resized to fit the new groundwater table, by means of changing the cell
% sizes of each UZ-model, while preserving the pressure profile. No iterations
% are  performed here, which makes the model rather inexact and approximate 
% (and by no means mass-conservative), but very very fast. Hence, USE WITH CARE! 
%
%---------------------------------------------------------------------
% Input:  UZ - structure containing all UZ-model information
%         GRW - structure containing all GRW-model information
%         COUP- structure containing all coupling-related info
%         hUZorg - cell with initial UZ heads (one for each UZ-model, all
%               ranging from bedrock to landsurface
%         hGRW - matrix with initial groundwater heads
%         sDir - directory to save the results in
%         reStartFrom - restart position indentifier
%
%---------------------------------------------------------------------
% LATEST EDITION
% Stripped version created for the publishing
%
%---------------------------------------------------------------------

function runCoupled_V3p2(UZ,GRW,COUP,hUZorg,hGRW,sDir,reStartFrom)
tic

addpath grwCPU
addpath UZ-ISU1Dv4

disp('Welcome to the Coupled UZ-GRW model')
disp('This is coupled version 3.2.2 (D. Erdal, Uni Tübingen, 2016/2017)')
disp(' ')
disp(['Resutls will be save in: ' sDir ])
disp(' ')
tTot=GRW.grwms.tTot;
dt=GRW.grwms.dt;


% Sort out a few technical things
UZ.delzCS=UZ.geom.delzCS;
UZ.delzNN=UZ.geom.delzNN;
UZ.solver.doEnKF=true;
GRW.nrCell=numel(GRW.grwms.K);
GRW.solver.doPrint=false;


% Create saving directory and save the initial setup
mkdir(sDir)
save([sDir '/Temp_START'])

% Small sanity check before starting
% Number of cells && Timing
if numel(GRW.grwms.K)<size(UZ.kin,2) || UZ.time.dt > GRW.grwms.dt
    warning('THE SETUP IS STRANGE --> CHECK SETUP!!!')
end

%% Coupled model initialization

% number of unsaturated zone columns
nrUZ=size(hUZorg,2);
% UZ-Index
ixx = false(size(hUZorg));
% Cell sizes for the unsaturated zone models
delzCS=cell(nrUZ,1);
% Water content before resizing
wcOLD=cell(nrUZ,1);
% Pressure head in UZ
hUZ=cell(nrUZ,1);
% Current water content in UZ
wcUZ=cell(nrUZ,1);

% Counter for number of cells (nodes) in each UZ column
nrCells=zeros(nrUZ,1);

% Loop through each of the UZ columns and re-construct such that each
% column ranges from the groudnwater table to the surface
% (originally each column goes until bedrock)
for i=1:nrUZ
    
    % find the interface between the current UZ and its corresponding GRW
    tmp=find(hUZorg(:,i)<0);
    ix=tmp(1);
    
    % find the exact placement of the water table
    dzGRW=interp1(hUZorg(ix-1:ix,i),[0,UZ.delzNN(ix)],0)+UZ.delzCS(ix-1)/2;
    
    if abs(dzGRW-UZ.delzCS(ix-1)) < 2*eps
        % h=0 is at a boundary interface
        ixx(ix:end,i)=true;
        nrCells(i)=sum(ixx(:,i));
        delzCS{i}=UZ.delzCS(ixx(:,i));
        hUZ{i}=hUZorg(ixx(:,i),i);
        
        % ELSE: assign a smaller bottom cell such that h=0 is at the interface
    elseif dzGRW > UZ.delzCS(ix-1)
        % the groundwater table is in the ix cell
        % --> decrease that cells size
        ixx(ix:end,i)=true;
        nrCells(i)=sum(ixx(:,i));
        delzCS{i}=UZ.delzCS(ixx(:,i));
        hUZ{i}=hUZorg(ixx(:,i),i);
        
        % tune the cell size and bottom boundray
        delzCS{i}(1) = (UZ.delzCS(ix-1)+UZ.delzCS(ix)) - dzGRW;
        hUZ{i}(1)=interp1([0,delzCS{i}(1)+delzCS{i}(2)/2],[0,hUZorg(ix+1,i)],delzCS{i}(1)/2);
        
        
    elseif dzGRW < UZ.delzCS(ix-1)
        % groundwater tabel is in ix-1
        ixx(ix-1:end,i)=true;
        nrCells(i)=sum(ixx(:,i));
        delzCS{i}=UZ.delzCS(ixx(:,i));
        hUZ{i}=hUZorg(ixx(:,i),i);
        
        % tune the cell size and bottom boundray
        delzCS{i}(1)= delzCS{i}(1)-(delzCS{i}(1)/2 + hUZorg(ix-1,i));
        hUZ{i}(1)=interp1([0,delzCS{i}(1)+delzCS{i}(2)/2],[0,hUZorg(ix,i)],delzCS{i}(1)/2);
        
    else
        error('This should not happen.....')
    end
    
end

% Height between bedrock and surface
heightTot=UZ.geom.heighcell;
% Orignal vertical UZ grid
zORG=zeros(length(UZ.delzCS),1);
zORG(1)=UZ.delzCS(1)/2;
for j=2:length(UZ.delzCS)
    zORG(j)=(UZ.delzCS(j-1)+UZ.delzCS(j))/2;
end
zORG=cumsum(zORG(1:end));


% Create a copy of the UZ structure that does not contain the original
% parameter matrices (kin, aplpha, etc)
UZf.param=UZ.param;
UZf.boundary=UZ.boundary;
UZf.geom=UZ.geom;
UZf.time=UZ.time;
UZf.solver=UZ.solver;
% Extract the UZ parameter matrices
KINX=UZ.kin;
ALPHAX=UZ.alpha;
NX=UZ.n;
SRESX=UZ.Sres;
SSATX=UZ.Ssat;
POROX=UZ.poro;
AKRX=UZ.aKR;

% Check if we are doing a restart or starting from scratch
if ~isempty(reStartFrom)
    strt=reStartFrom+1;
    % Load previously generated data
    load([sDir '/Temp_' num2str(reStartFrom)],'hGRW','hUZ','hGRWold','delzCS','R','R_rest','R_UZ','wcOLD');
else
    strt=1;
end

% Cumulative rechare (for error checking)
R_cum=zeros(nrUZ,1);
% Initilize the temporary parameter saver
ptmp=cell(nrUZ,1);

%% Start the General Outer Loop of Coupling Time Steps
% (here: the same as the groundwater time steps)

prIndex=max(1,round(round(tTot/dt)/100)); % print-index
prcount=0; % print-counter
disp('====Starting calculations========')

for i=strt:round(tTot/dt)
    % Initilize the recharge vector from UZ to GRW
    R_UZ=zeros(nrUZ,1);
    
    % Print
    if mod(i,prIndex)==0
        fprintf([num2str(round(100*100*i/(tTot/dt))/100) ' '])
        if prcount == 9; fprintf(1,'\n'); prcount=0;
        else, prcount=prcount+1;
        end
    end
    
    %% 1) Do the UZ
    UZf.time.tStop=i*dt;
    UZf.time.tStart=(i-1)*dt;
    
    if isfield(UZ.boundary,'fluxTopMTRX')
        % cut out only the useful flux-times (reuces paralellization overhead for long
        % simulations)
        if1=find(min(abs(UZ.boundary.changeTimes-UZf.time.tStart))==abs(UZ.boundary.changeTimes-UZf.time.tStart));
        if2=find(min(abs(UZ.boundary.changeTimes-UZf.time.tStop))==abs(UZ.boundary.changeTimes-UZf.time.tStop));
        
        UZf.boundary.fluxTopMTRX=UZ.boundary.fluxTopMTRX(:,max(1,if1-5):min(size(UZ.boundary.fluxTopMTRX,2),if2+5));
        UZf.boundary.changeTimes=UZ.boundary.changeTimes(:,max(1,if1-5):min(size(UZ.boundary.fluxTopMTRX,2),if2+5));
    end
    
    % Loop through each of the unsaturated columns
    % --> PARFOR
    for j=1:nrUZ
        
        % Get the base model structures
        param=UZf.param;
        boundary=UZf.boundary;
        geom=UZf.geom;
        time=UZf.time;
        solver=UZf.solver;
        
        % Get number of active cells
        geom.nz=nrCells(j);%
        
        % Get/Compute the z vectors
        geom.delzCS=delzCS{j};
        z=zeros(length(delzCS{j})+1,1);
        z(1)=delzCS{j}(1)/2;
        for jj=2:length(delzCS{j})
            z(jj)=(delzCS{j}(jj-1)+delzCS{j}(jj))/2;
        end
        z(end)=delzCS{j}(end)/2;
        geom.delzNN=z;
        z=heightTot-sum(z)+cumsum(z(1:end-1));
        
        % Total height of the current cell
        geom.heighcell=sum(geom.delzCS);
        
        % Get the right parameters by interpolation from the original
        % parameters
        param.kin = interp1(zORG,KINX(:,j),z,'linear','extrap');
        param.alpha = interp1(zORG,ALPHAX(:,j),z,'linear','extrap');
        param.n = interp1(zORG,NX(:,j),z,'linear','extrap');
        param.Ssat = interp1(zORG,SSATX(:,j),z,'linear','extrap');
        param.Sres = interp1(zORG,SRESX(:,j),z,'linear','extrap');
        param.poro = interp1(zORG,POROX(:,j),z,'linear','extrap');
        param.aKR = interp1(zORG,AKRX(:,j),z,'linear','extrap');
        param.m=1-1./param.n;
        
        
        % Assign initial condition
        boundary.headinit=hUZ{j};%(ix,j);
        
        % Assign the flux bounary if not already done
        if isfield(boundary,'fluxTopMTRX')
            % find the zone number
            iflx=COUP.zoneMap(COUP.UZindex(j,1),COUP.UZindex(j,2));
            boundary.fluxTop=boundary.fluxTopMTRX(iflx,:);
        end
        
        % Run the UZ model
        % Ouput: head, recharge, water content
        [hUZ{j},R_UZ(j,1), wcUZ{j}]=solverStaticCOUP_app2(geom,param,solver,time,boundary);
        
        % ad-hoc rescaling of R to account for potential differences between the
        % porostiy in the UZ and the specific yield in the GRW
        R_UZ(j,1)= R_UZ(j,1)/param.poro(1)*GRW.grwms.poro(COUP.UZindex(j,1),COUP.UZindex(j,2));
        
        
        if any(isinf(hUZ{j}))
            error(['UZ solver crashed for zone number ' num2str(j)])
        end
        
        % and get the water content in the system
        param.sDiff=param.Ssat-param.Sres;
        wcOLD{j}=calcSat(hUZ{j},param).*param.poro;
        % save the current parameter setups
        ptmp{j}=param;
    end
    
    %% Regionalize the UZ recharge based on the zone map
    if numel(GRW.grwms.K)~=size(UZ.kin,2)
        if max(max(COUP.zoneMapIndivid))~=size(UZ.kin,2)
            warning('OBS: NOT IMPLEMENTED THIS WAY!!!')
        end
        
        % interpolate  according to zoneMapIndivid
        R=zeros(GRW.grwms.ny,GRW.grwms.nx);
        
        for iz=1:max(max(COUP.zoneMapIndivid))
            % one UZ column per individual land-use
            R(COUP.zoneMapIndivid==iz)=mean(R_UZ(iz));
        end
        
    else
        % Each cell in the GRW has its own UZ column, no interpolation
        % needed
        R=R_UZ;
    end
    
    % cummulative recharge (for error search etc.)
    R_cum=R_cum+R_UZ;
    
    
    %% 2) Run the GRW-model
    
    GRW.grwms.tTot=i*dt; % Stopping time
    GRW.grwms.tStart=(i-1)*dt; % Starting time
    GRW.boundary.Q=-reshape(R/dt,GRW.grwms.ny,GRW.grwms.nx); % Recharge in correct shape and units
    
    GRW.grwms.hInit=reshape(hGRW,GRW.grwms.ny,GRW.grwms.nx); % Initial head
    hGRWold=reshape(hGRW,GRW.grwms.ny,GRW.grwms.nx);    % for later use
    
    % Run the model. Output: head, structure info
    [hGRW,~]=runGRWunconf(GRW.solver,GRW.grwms,GRW.boundary);
    
    if any(any(isinf(hGRW)))
        error('GRW-model crashed')
    end
    
    % Reshape
    hGRW=reshape(hGRW,GRW.grwms.ny,GRW.grwms.nx);
    
    %% 3) Post-process the intput and redefine the grid
    for j=1:nrUZ
        
        % Get the difference in groundwater head for this timestep
        if numel(GRW.grwms.K)~=size(UZ.kin,2)
            dGRW=hGRW(COUP.UZindex(j,1),COUP.UZindex(j,2))-hGRWold(COUP.UZindex(j,1),COUP.UZindex(j,2));
        else
            dGRW=hGRW(j)-hGRWold(j);
        end
        
        % if change is almost zero, let it be 0
        if abs(dGRW) < 3*eps
            dGRW=0;
        end
        
        % If the extent is too small (<2cm), then use this extent
        if GRW.grwms.z0(COUP.UZindex(j,1),COUP.UZindex(j,2)) + UZ.geom.heighcell - hGRW(COUP.UZindex(j,1),COUP.UZindex(j,2)) < 0.02
            delzCS{j}=delzCS{j}/sum(delzCS{j})*0.02;
            dGRW=0;
        end
        
        % adjust the grid so that the number of cells remains but their
        % sizes are adjusted to account for the new UZ zone size
        a=delzCS{j}/sum(delzCS{j})*dGRW;
        % new grid
        delzCS{j}=delzCS{j}-a;
        
        if any(delzCS{j}<0)
            error('!!! NEGATIVE CELL SIZES !!!!')
        end
    end
    
    
    
    %% 4) Check if a Saving is Due
    if mod(i,COUP.saveEach)==0
        save([sDir '/Temp_' num2str(i)],'hGRW','hUZ','hGRWold','delzCS','R','R_rest','R_UZ','wcOLD','R_cum','X','wcUZ')
        R_cum=zeros(nrUZ,1);
    end
    
end

%% Finish the program
% Save a dummy variable so that the start-file is listed at the top
dummy=nan;
pause(1)
save([sDir '/Temp_START'],'dummy','-append')

if prcount ~= 0; fprintf(1,'\n'); end
disp('====Calculations Finished========')

toc
