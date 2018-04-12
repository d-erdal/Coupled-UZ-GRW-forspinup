%% PROGRAM INFORMATION
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date: August/September 2014,
% GPU version March/April 2015
%
% The program solves the nonlinear equations for a transient unconfined (phreatic)
% aquifer in horizontal 2D: n*(dh/dt) - div(K*(h-z)*grad(h)) = Q
% where n is porosity, K saturated conductivity, ground elevation, Q
% recharge and h groundwater level over NN.
%
% Discretization: Finite volume in space with fully implicit Euler in time
% Semi-discretized (dx=dy, sum of neighbours j to the cell i):
% n*(hi-hi,old)/dt - (1/(2dx²))*sum(K*(hj²-hi²+2*z*(hi-hj)) - Qi = 0
%
% Output:   hStore - head values at selected positions
%           extra - additional inforation including the indexing for use
%           with a non-rectangular domain


% LATEST EDITION
% Stripped version created for the publishing
%


%%
function [hStore,extra]=runGRWunconf(solver,grwms,boundary)

%% Initialize the system
extra=[]; ticTot=tic; cTot=0; hStore=[];


%% Welcome text

if solver.doPrint
    str=['\n\n\n====' date '========='];
    str=str(1:28);
    str=[str, '\nWelcome to the unconfined aquifer solver! \nThis is version CPU 2.1, 2015. \n======================'...
        '\nProblem run time: ' num2str(grwms.tTot)...
        ' s (~' num2str(round(grwms.tTot/3600)) 'h, ~' num2str(round(grwms.tTot/3600/24)) 'd) with time step ' ...
        num2str(grwms.dt) ' s (~' num2str(round(grwms.dt/3600)) ' h, ~' num2str(round(grwms.dt/3600/24)) ' d).\n'];
    fprintf(str)   
   str='Solving with Newtons method and ';   
    
    if solver.how == 1
        str = [str 'Matlab backslash operator'];
    elseif solver.how == 2
        str = [str 'Biconjugate gradients stabilized method (bicgstab)'];
    elseif solver.how == 3    
        str = [str 'Jacobis method on CPU'];
    end
    
    if solver.doExplicitPredcond
        str = [str '\nPreconditioned by Explicit Euler'];
    else
        str = [str '\nNo preconditioning'];
    end
    
    if solver.doLineSearch
        str = [str ' and using LineSearch'];
    end
    
    str=[str,'\n======================\n'];
    fprintf(str)
end
%% Sanity checks
% Check if there are problems with the setup before starting
a=grwms.K==0 & boundary.mask~=0;
b=grwms.poro==0 & boundary.mask~=0;


if any(any(b)) || any(any(a))
    if any(any(b)) && grwms.dt==1 && grwms.tTot==1
        if solver.doPrint
            disp('Running for steady state!')
        end
        solver.doSS=true;
    else
        warning('Some parameters are 0 --> Check problem setup!')
    end
end
clear a b

% For compability with older versions:
% try    boundary.QtempoScale;
% catch; boundary.QtempoScale=1;
% end

if boundary.doRiver
if ~isfield(boundary,'boundary.hRivertempoScale')
boundary.hRivertempoScale=1;
end
end

if ~isfield(solver,'doSS')
solver.doSS=false;
end

if ~isfield(grwms.wells,'newAppr')
    if isfield(grwms.wells,'changeW') 
        grwms.wells.newAppr=true;
    else grwms.wells.newAppr=false;
    end
end

try
    if solver.saveToDisk
        mkdir(solver.dirName)
    end
catch
    solver.saveToDisk=false;
end

%% Practical variables
nx=size(boundary.mask,2);   % number of cells in x
ny=size(boundary.mask,1);   % number of cells in y
dy=grwms.dy;                % cell sizes
dx=zeros()+grwms.dx;                % cell sizes
dt=grwms.dt;                % time step

if dx~=dy
    % this has to do with the mdxdy parameter in the matrix creation
    
%OPEN QUESTION SINCE PREVIOUS: IS THIS STILL TRUE? I GUESS NOT.....
% --> YES: VERY MUCH SO (WHY??) (08.09.2016)
    warning('Diriclet boundaries may be wrong if dx and dy are strongly different')
end
%% Sort out indexes:
% Welcome to an indexing nightmare; this takes into account also
% non-rectangular systems.

% Total index (index of each active cell)
tmp=reshape(1:nx*ny,ny,nx);
tmp(boundary.mask==0)=nan;
tmp=reshape(tmp,nx*ny,1);
tmp(isnan(tmp))=[];
indexTot=tmp;

% Number of active cells:
nrCells=length(indexTot);

%===Neighbour index (top, bottom, left, right)===

tmp=nan(ny,nx);
tmp(boundary.mask==1)=1:nrCells;
tmp2=nan(ny+2,nx+2);
tmp2(2:end-1,2:end-1)=tmp;

% Actual cell index
ime(1:nrCells,1)=1:nrCells;

% Right hand neighbour
tmp=tmp2(2:end-1,3:end);
tmp(boundary.mask==0)=nan;
ir=tmp(indexTot);       % the neighbours index
irDo=true(size(ir));    % index as boolean (used to create the Jacobian)
irDo(isnan(ir))=false;
ir(isnan(ir))=1;

% Left hand neighbour
tmp=tmp2(2:end-1,1:end-2);
tmp(boundary.mask==0)=nan;
il=tmp(indexTot);
ilDo=true(size(il));
ilDo(isnan(il))=false;
il(isnan(il))=1;

% Top neighbour
tmp=tmp2(3:end,2:end-1);
tmp(boundary.mask==0)=nan;
it=tmp(indexTot);
itDo=true(size(it));
itDo(isnan(it))=false;
it(isnan(it))=1;

% Bottom neighbour
tmp=tmp2(1:end-2,2:end-1);
tmp(boundary.mask==0)=nan;
ib=tmp(indexTot);
ibDo=true(size(ib));
ibDo(isnan(ib))=false;
ib(isnan(ib))=1;

% What goes where in the Jacobian matrix
iA1=([ime;ime(itDo);ime(ibDo);ime(ilDo);ime(irDo)]);
iA2=([ime;it(itDo);ib(ibDo);il(ilDo);ir(irDo)]);

% the index of the perimiter cells (boundary cells)
iPerim=boundary.perimI(indexTot);
extra.indexTot=indexTot; % save the index for post-processing
clear tmp tmp2

%% Make sure that a boundary is only on a cell with a boundary!
% (hence: delete from the boundary index cells that are on inner edges on a boundary)
maskL=zeros(ny+2,nx+2);
maskL(2:end-1,2:end-1)=boundary.mask;
b1=diff(maskL(:,2:end-1),1,1)~=0;
b2=diff(maskL(2:end-1,:),1,2)~=0;
clear maskL

% find out how many boundaries each cell holds
% (for the GPU implementation this is not used)
perim=(b1(1:end-1,:)+b1(2:end,:)+b2(:,1:end-1)+b2(:,2:end)).*boundary.mask;
iPerim(perim(indexTot)==0)=0;
clear perim b1 b2

%% Pre-process the river boundary
% Asumption: river is on a finer (or same) grid as the problem AND the two
% grids are just different scales of each other (dx/riverDX=integer>=0)
if boundary.doRiver
    r=grwms.dx/boundary.riverDX;
    if r==1
        Ar=boundary.riverDX*boundary.riverDY*boundary.riverMask;
    elseif isfield(boundary,'riverAr') && all(size(boundary.riverAr)==[ny,nx])
        Ar=boundary.riverAr;
    else        
        Ar=zeros(size(boundary.mask));
        % Loop through each actual grid cell
        % (could be speeded up on the gpu.....)
        for i=1:ny
            for j=1:nx
                Ar(i,j)=sum(sum(boundary.riverMask(1+(i-1)*r:i*r,1+(j-1)*r:j*r)));
            end
        end
        Ar=boundary.riverDX.*boundary.riverDY.*Ar;
    end
    Ar(boundary.mask==0)=0;
 %   Ar(boundary.perimI~=0)=0;
    
    KRiver=(boundary.Kriver./boundary.dRiver.*Ar./(grwms.dx*grwms.dy));
    Ar(Ar~=0)=1;
    boundary.riverMaskDomain=logical(Ar);
    clear Ar
    
    KRiverI=KRiver(indexTot);
else
    KRiverI=0;
end


%% Wells: part 1
% Get the index of the cells with the wells

if grwms.doWells
    wXY=zeros(ny,nx);
    for i=1:grwms.wells.nr
        wXY(ceil(grwms.wells.y(i)/dy),ceil(grwms.wells.x(i)/dx))=i;
    end
else
    wXY=0;
end

%% Conductivity (K) and bedrock elevation (Z)
% Get the values at the interfaces between the cells
%

grwms.K(boundary.mask==0)=nan;
K=(nan(ny+2,nx+2));
K(2:end-1,2:end-1)=grwms.K;
% Call the GPU function for averaging the Ks
dxF=zeros(ny,nx)+dx; dyF=zeros(ny,nx)+dy; 
[ki,kt,kb,kl,kr]=elementK_CPU(K(2:end-1,2:end-1),K(3:end,2:end-1),K(1:end-2,2:end-1),K(2:end-1,1:end-2),K(2:end-1,3:end),dxF,dyF);
clear K

ki=ki(indexTot);    % Conductivity for the potential Diricelt boundary
kt=kt(indexTot);    % Conductivity at top interface
kb=kb(indexTot);    % Conductivity at bottom interface
kl=kl(indexTot);    % Conductivity at left interface
kr=kr(indexTot);    % Conductivity at right interface

% Same as for K but for Z
Z=(nan(ny+2,nx+2));
Z(2:end-1,2:end-1)=grwms.z0;
[zi,zt,zb,zl,zr]=elementZ_CPU(Z(2:end-1,2:end-1),Z(3:end,2:end-1),Z(1:end-2,2:end-1),Z(2:end-1,1:end-2),Z(2:end-1,3:end));
clear Z

zi=zi(indexTot);
zt=zt(indexTot);
zb=zb(indexTot);
zl=zl(indexTot);
zr=zr(indexTot);

% Porosity on the GPU
poro=(grwms.poro(indexTot));

%% Initilize boundaries, recharge and rivers


% ---------------- BOUNDARIES ------------------------------------------
boundNeu=zeros(nrCells,1);
boundDir=nan(nrCells,1);

% Initialize the Diriclet boundaries
for j=1:length(boundary.perimIdir)
    boundDir(iPerim==boundary.perimIdir(j))=boundary.boundDirValues(j,1);
end
% -----------------------------------------------------------------------
% Initialze the Neumann boundaries
for j=1:length(boundary.perimIneu)
    boundNeu(iPerim==boundary.perimIneu(j))=boundary.boundNeuValues(j,1);
end
% -----------------------------------------------------------------------

% ---------------- RECHARGE----------------------------------------------
if isfield(boundary,'Qvector') && boundary.doQvector
    Quse=zeros(size(boundary.Qmap));
    for imap=1:size(boundary.Qvector,1)
        Quse(boundary.Qmap==imap)=boundary.Qvector(imap,1);
    end    
else
    Quse=(boundary.Q(:,:,1)*boundary.QtempoScale(1));
end
Qorg=boundary.Q(:,:,1);
Qorg=Qorg(indexTot);
% -----------------------------------------------------------------------


% ---------------- RIVERS----------------------------------------------
if boundary.doRiver
    % Initialize the river stage value (for placement see above)
    if any(size(boundary.hRiver)==1)
        % inpuit is a scalar --> same stage at full river
        hRuse=zeros(ny,nx) + boundary.hRiver(1)*boundary.hRivertempoScale(1) + grwms.lsurf;
        if length(boundary.hRiver)==1
            % Potential temporal changes goes in hRivertempoScale
            hRwht=1;
        else
            % Temporal changes are in boundary.hRiver
            hRwht=2;
        end
    else
        % input is distributed --> different cells have differnt values
        hRuse=(boundary.hRiver(:,:,1)*boundary.hRivertempoScale(1) + grwms.lsurf);
        if size(boundary.hRiver,3)==1
            % Potential temporal changes goes in hRivertempoScale
            hRwht=1;
        else
            % Temporal changes are in boundary.hRiver
            hRwht=3;
        end
        
    end
    hRuse(KRiver==0)=0;
    
    if hRwht==1
        % If hRivertempoScale is to be used; get the original value
        hRorg=(hRuse-grwms.lsurf)./boundary.hRivertempoScale(1);
    end
    
else
    hRuse=(0);
end
% -----------------------------------------------------------------------

%% Setup for Transient Solution

% Create the storage vector for the pressure heads
if ~solver.saveToDisk
    hStore=zeros(nrCells,ceil((grwms.tTot-grwms.tStart)/grwms.dtStore));
end
extra.tStore=zeros(1,ceil((grwms.tTot-grwms.tStart)/grwms.dtStore)); % (actual times for storage)

% Old head = inital head
if numel(grwms.hInit)==nx*ny
	hOld=(grwms.hInit(indexTot));
else
	hOld=grwms.hInit;
end

if solver.doPrint
    disp('Now starting the calculations (% done):')
end

% initialize solver parameters, counters etc
zahler2=0; zahler4=0; ti=grwms.tStart; stIndex=round(ti/grwms.dtStore)+1; stOrg=stIndex-1;
forceQuit=false; prcount=-1; tscount=0; prIndex=round(grwms.tTot/dt/100);
if prIndex==0; prIndex=1; end

%%
tct=tic;
extra.QTS=[];
% Start the time loop
while ti<grwms.tTot
    
    
    % Print current status
    if solver.doPrint && mod(tscount,prIndex)==0
        fprintf([num2str(round(100*100*ti/grwms.tTot)/100) ' '])
        if prcount == 9; fprintf(1,'\n'); prcount=0;
        else prcount=prcount+1; end
    end
    
    % Update the Diriclet boundaries
    if boundary.chngDir
        % find the right index
        tTMP=ti-boundary.changeTdir;
        tTMP(tTMP<0)=NaN;
        
        for j=1:length(boundary.perimIdir)
            boundDir(iPerim==boundary.perimIdir(j))=boundary.boundDirValues(j,min(tTMP)==tTMP);
        end
    end
    
    % Update the Neumann boundaries
    if boundary.chngNeu
        % find the right index
        tTMP=ti-boundary.changeTneu;
        tTMP(tTMP<0)=NaN;
        
        for j=1:length(boundary.perimIneu)
            boundNeu(iPerim==boundary.perimIneu(j))=boundary.boundNeuValues(j,min(tTMP)==tTMP);
        end
    end
    
    % For recharge, wells and rivers: fill out the Q matrix that goes on
    % the right hand side
    
    % Update the Reacharge
    if boundary.chngQ
        tTMP=ti-boundary.changeTQ;
        tTMP(tTMP<0)=NaN;
        if size(boundary.Q,3)>1
            % Q is a nx x ny x changeTimes matrix
            Q=(full(Qvec(indexTot,min(tTMP)==tTMP)));
            % 08.07.16 --> does this work in the current version?
        elseif isfield(boundary,'Qvector') && boundary.doQvector
            % implemented 08.07.16
            Q=zeros(size(boundary.Qmap));
            for imap=1:size(boundary.Qvector,1)
                Q(boundary.Qmap==imap)=boundary.Qvector(imap,min(tTMP)==tTMP);
            end
            Q=Q(indexTot);
        else
            % Q is a sacling field scaled by QtempScale
            Q=Qorg*boundary.QtempoScale(min(tTMP)==tTMP);
            % extra.QTS=[extra.QTS;boundary.QtempoScale(min(tTMP)==tTMP)];
        end
    else
        Q=Quse(indexTot);
    end
    
    % Update the Wells
    if grwms.doWells
        if grwms.wells.newAppr
            % New approach!
            tTMP=ti-grwms.wells.changeW;
            tTMP(tTMP<0)=NaN;
            W=zeros(ny,nx);
            for iw=1:grwms.wells.nr
                W(wXY==iw)=grwms.wells.Q(min(tTMP)==tTMP,iw);
            end
            
        else
            % if not in the start-file, use old approach
            
            W=zeros(ny,nx);
            % Check if any wells are active
            activeW=ti-grwms.wells.tStart >= 0 & ti-grwms.wells.tStop < 0;
            if any(activeW)
                for iw=1:grwms.wells.nr
                    if activeW(iw)
                        W(wXY==iw)=grwms.wells.Q(iw);
                        % grwms.wells.Q(iw)
                    end
                end
            end
        end
        Q=Q+W(indexTot);
    end
    
    % Update the River boundary
    if boundary.doRiver 
        if boundary.chngRiver
            tTMP=ti-boundary.changeTRiver;
            tTMP(tTMP<0)=NaN;
            
            if hRwht==1
                hRuse=hRorg*boundary.hRivertempoScale(min(tTMP)==tTMP) + grwms.lsurf;
            elseif hRwht==2
                hRuse=zeros(ny,nx) + boundary.hRiver(min(tTMP)==tTMP) + grwms.lsurf;
                hRuse(KRiver==0)=0;
            elseif hRwht==3
                hRuse=(boundary.hRiver(:,:,min(tTMP)==tTMP)+ grwms.lsurf);
                hRuse(KRiver==0)=0;
            end
        end
        Q=Q-KRiverI.*hRuse(indexTot);
    end
    
    
    %% PRECONDITION
    if solver.doExplicitPredcond
        % Solve the time step with explicit Euler and use the soltion
        % as initial quess for the nonlinear iterations.
        hOldE=hOld;
        inflow=(-Q + boundNeu);
        tx=0;
        c=1;
        while tx<grwms.dt
            
            [flowB,dtMin]=elementEulerExpl_CPU(hOldE(ime),hOldE(it),hOldE(ib),hOldE(il),hOldE(ir),kt,kb,kl,kr,zt,zb,zl,zr,boundDir,ki,zi,inflow,poro,KRiverI,solver.EulerMax);
            dtX=min(min(dtMin),grwms.dt-tx);
            hOldE=hOldE+dtX.*flowB;
            tx=tx+dtX;
            c=c+1;
            
            if c>10000
                hOldE=hOld;
                break
            end
            
        end
        hItOld=(hOldE);
     %   clear hOldE flowB dtMin
    else
        hItOld=(hOld);
    end
    
    dtF=zeros(nrCells,1)+dt;
    %% START THE SOLVER
    zahler3=0;   zMax=100; if solver.doSS; zMax=zMax*2; end
    
    while 1
        
        dh=inf; zahler1=0; normf1=inf;
        
        %% Start of the nonlinear iterations
        hOld=(hOld);   % (assure it is on the GPU)
        
        while max(abs(dh))>1e-8 && zahler1<zMax && (normf1>1e-7 || max(abs(dh))>1e-5)
           
            if solver.doSS && length(hOld)>1e6
                % we have a slow problem, display convergance features
                disp(max(abs(dh)))
            end
            
            inflow=(-Q + boundNeu);
            
            
          % Compute the entries of the Jacobian             
            [ji,jt,jb,jl,jr,f0]=elementJacobian_CPU(hItOld(ime),hItOld(it),hItOld(ib),hItOld(il),hItOld(ir),...
                kt,kb,kl,kr,zt,zb,zl,zr,poro,boundDir,ki,zi,hOld,inflow,dtF,KRiverI);
%             wait(gd) % Include if using the profiler, otherwise the times are looking odd
            
            if solver.how == 1
%================= use the MATLAB backslash operator=======================
                
                % Construct the Jacobian
                J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));(jr(irDo))]);
                % Solve the system
                dh=J\-(f0);                
                
            elseif solver.how == 2
%================= use the iterative BICGSTAB METHOD=======================

                % Construct the Jacobian
                J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));(jr(irDo))]);
                % Incomlete LU-factorization for preconditioning
                [L,U]=ilu(J);
                % Solve the system
                [dh,out]=bicgstabl(J,-(f0),solver.cErr,2000,L,U);
                
                % Looking like the system did not converge
                if out~=0
                    % check if we actually have a problem
                    if max(abs(J*dh+f0))>1e-9
                        % We do: solve with \
                        disp('No covergence, solving using mldivide')
                        dh=J\-(f0);
                    end
                end
                
            elseif solver.how == 3
%================= use the iterative JACOBI METHOD=========================                
                
               % Initial guess (0 is here often best)
                dx=zeros(nrCells,1);
                
              c=1;
              
              
              % Current setup: calculate the error in the same way as the
              % BICGSTAB.... Rather costly but gives comparable results
              
              % Norm of the right hand side
              normf0=norm(f0,2);
              while 1
                  % Compute the new dx
                  [dxnew]=solveEQ_CPU(ji,jt,jb,jl,jr,dx(it),dx(ib),dx(il),dx(ir),-f0);
                  %                     wait(gd);
                  
                  
                  
                  % Calculate the norm (sum of a GPU-array is EXPENSIVE!)
                  eX=forNorm_CPU(ji,jt,jb,jl,jr,dxnew,dxnew(it),dxnew(ib),dxnew(il),dxnew(ir),f0);
                  err2=sqrt(sum(eX))/normf0; % norm(f0-J*dxnew)/norm(f0)
                  
                  
                  %  err=((max(abs(dxnew-dx)))); % Alternative error that is cheeper to compute but less accurate
                  if err2 < solver.cErr
                      % Convergance: save iteration number and exit
                      cTot=cTot+c;
                      break
                      
                  elseif c>solver.maxJacobi
                      % If no convergance, use \
                      disp('No covergence, solving using mldivide')
                      J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));(jr(irDo))]);
                      dxnew=J\-(f0);
                      break
                  end
                  % Iterate on.....
                  c=c+1;
                  dx=dxnew;
              end
              clear dx
              dh=dxnew;
              clear dxnew
              % Update counter
              zahler4=zahler4+c;
            
            
             elseif solver.how == 4
%================= MIXTURE =======================
            if zahler1<=1
                % do matlab
                 % Construct the Jacobian
                J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));(jr(irDo))]);
                % Solve the system
                dh=J\-(f0);                
            else
%                 % do bicgstab
%                 
                % Construct the Jacobian
                J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));(jr(irDo))]);
                % Incomlete LU-factorization for preconditioning
                [L,U]=ilu(J);
                % Solve the system
                [dh,out]=bicgstab(J,-(f0),solver.cErr,200,L,U);
                
                % Looking like the system did not converge
                if out~=0
                    % check if we actually have a problem
                    if max(abs(J*dh+f0))>1e-9
                        % We do: solve with \
                        disp('No covergence, solving using mldivide')
                        dh=J\-(f0);
                    end
                end


            end
            end
         %   clear ji jt jb jl jr
%==========================================================================

            %% DO LINESEARCH            

            alpha=1;
            if solver.doLineSearch
                
                % original norm
                if solver.how~=3                   
                   normf0=norm(f0,2);
                end
                
                % alpha=1 norm
                hX=hItOld+dh;
                f1_1=elementf0_CPU(hX(ime),hX(it),hX(ib),hX(il),hX(ir),...
                    kt,kb,kl,kr,zt,zb,zl,zr,poro,boundDir,ki,zi,hOld,inflow,dtF,KRiverI);
                normf1_1=norm(f1_1,2);
                normf1=normf1_1;
                
                % Iterate until the new soltuion has a smaller norm than
                % the old one
                while normf1 >= normf0 && alpha > 2*1/1000 && normf1>1e-6
                    alpha=alpha*0.5;
                    hX=hItOld+alpha*dh;
                    f1=elementf0_CPU(hX(ime),hX(it),hX(ib),hX(il),hX(ir),...
                        kt,kb,kl,kr,zt,zb,zl,zr,poro,boundDir,ki,zi,hOld,inflow,dtF,KRiverI);
                    normf1=norm(f1,2);
                end
                
                if alpha <= 2*1/1000 %(max(abs(dh))*alpha) <= 0.01%
                    % hence, if nothing good was found
                    % take the largest or smallest step
                    if round(normf1/normf0) >= round(normf1_1/normf0)
                        normf1=normf1_1;
                        alpha=1;
                    end
                end
                
            else
                % Just calculate the norm of the new soltion
                % (if needed for judging convergance)
%                 hX=hItOld+alpha*dh;
%                 f1=elementf0_CPU(hX(ime),hX(it),hX(ib),hX(il),hX(ir),...
%                     kt,kb,kl,kr,zt,zb,zl,zr,poro,boundDir,ki,zi,hOld,inflow,dtF,KRiverI);
%                 normf1=norm(f1,2);
                normf1=norm(f0,2);
                
            end
        %    clear hX f1 f1_1 f0 inflow
            
          %  hX=hItOld+alpha*dh;
            
           % TMP TMP TMP-----
%             alpha2=((0.01+zi)-hItOld)./(alpha*dh);
%             alpha2(alpha2>1)=1;            
%             disp(min(alpha2))
%             hNew=hItOld+max(min(alpha2),0.001)*alpha*dh;
            %-----------------
            
            % Update and go on
           hNew=hItOld+alpha*dh;
            
            % TMP TMP TMP------
%             hNew(hNew<zi)=inf;%zi(hNew<zi)+0.01;
             %-----------------
             
             
            zahler1=zahler1+1;
            hItOld=hNew;                    
            
        end
        
        
        
        %% Check the solution        
        
        % Update the counter
        zahler2=zahler2+zahler1;

        % Check if we are ok or if exited for wrong reasons
        if zahler1==zMax || any(isnan(hNew))
            
            if solver.doSS
                warning('--------->  No convergence reached')
                forceQuit=true;
                break
            end
            
            % No convergence: reduce time step and retry
            dt=dt/2;
            if solver.doPrint
                fprintf(['\nConvergence problem; new dt: ' num2str(dt) ' at time ' num2str(ti) ' ( ' num2str(zahler3) ' )\n'])
            end
            zahler3=zahler3+1; prcount=0;
            hItOld=(hOld);
            
        else
            % Update the time
            ti=ti+dt;
            break
        end
        
        if zahler3>10
            % after 10 times reduction of dt, quite
            forceQuit=true;
            break
        end
        
    end
       
    % if we broke out for real, then this is the end of the line: quit
    if forceQuit
        hStore=hStore+inf;
        break
    end
    
    
    %% STORE AND UPDATE
    
        
    if  grwms.dtStore*stIndex<=ti % check if this time step should be saved
        if solver.saveToDisk
            hSave=single(hNew);
            save([solver.dirName '/hNew_' num2str(stIndex)],'hSave')
        else      
            hStore(:,stIndex-stOrg)=(hNew);
        end
        extra.tStore(stIndex)=(ti);
        stIndex=stIndex+1;
    end
    
    % Old is new and on we go.....
    hOld=hNew;
    tscount=tscount+1;
 
end

% DONE!

%% Post Processing and printing
extra.ct=toc(tct);

if solver.doPrint
    fprintf(' 100!\n======================\n')
    if forceQuit
        error(['Severe convergence errors: The program is now exiting WITHOUT a complete solution. Exit time ', ...
            num2str(ti) ' out of ' num2str(grwms.tTot) '. Please reconsider solver settings and/or problem setup and restart'])
    else
        fprintf(['Problem solved using ' num2str(zahler2) ' number of nolinear iterations'])
        if solver.how==3
            fprintf(['\nand ' num2str(zahler4) ' Jacobi iterations.\n'])
        else
            fprintf('.\n')
        end
        
    end
    % toc
    toc(ticTot)
    fprintf('======================\n\n') 
end


extra.zNonLin=zahler2;
if zahler4~=0; extra.zJac=zahler4; end












