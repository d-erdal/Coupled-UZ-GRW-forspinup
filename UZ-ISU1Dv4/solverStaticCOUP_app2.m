% Part of a new version of the solution for the Richards equation
% Written Dec 2011 by Daniel Erdal
% Compacted version written May 2012

% Version 4.0
% Last Documented edition: 05.12.13 (DE)

% Note on this edition:
% The principle is the same as with version 3, but this version should work
% with both the EnKF and the Dream codes

function [headStore,R, wc_new]=solverStaticCOUP_app2(geom,param,solver,time,boundary) %geom,param,solver,time,boundary,tStart,tStop,headinter

R=0;

% geom.delzNN=zeros(geom.nz+1,1);
% geom.delzNN(1)=geom.delzCS(1)/2;
% for i=2:length(geom.delzCS)
%     geom.delzNN(i)=(geom.delzCS(i-1)+geom.delzCS(i))/2;
% end
% geom.delzNN(end)=geom.delzCS(end)/2;
% 
% geom.heighcell=sum(geom.delzCS);
% geom.nz=length(geom.delzCS);

%% -------------- Initialization --------------------------------------

if solver.doEnKF || solver.doDream
    doPrint=false;
else
    doPrint=solver.doPrint;
end

if doPrint
    fprintf(2,'ISU1Dv3: Newton-Raphson ')    %#ok<*PRTCAL>
    
    if time.doAdaptiveTime;                     fprintf(2,'with Adaptive time ');
    end
    if solver.doLineSearch;                 fprintf(2,'with Line Seach ')    
    end
    if boundary.doFreeDrain;                    fprintf(2,'with Free drain bottom\n')
    else                                        fprintf(2,'\n');
    end
end


% Calculate a few useful variables (just to not double their calculation later)
param.sDiff=param.Ssat-param.Sres;
param.nzKR=geom.nz+2;
boundary.KrelWeight=param.KrelWeight;
boundary.pWhat=param.what;

% Actual head variable: Initial condition
headinter=boundary.headinit;

% Old head variable: Actual head variable
head_old=headinter;

% for the free drain:
boundary.alphaB=param.alpha(1);
boundary.aKR=param.aKR(1);

% Initial (and old) saturation and water content
[sat,se]=calcSat(head_old,param);
wc_new=param.poro.*sat;

% set the grid
delzNN=geom.delzNN;
delzCS=geom.delzCS;

% initial time
tim=time.tStart;
zTot=0;

% time step counter
tscount = 1;

% print counter
prcount = 0;

% number of the output
feldcount = 1;

% storing variables
satStore=NaN;
headStore=NaN;
timestore=NaN;
tStore=[];

if ~solver.doEnKF
    if solver.doDream
        tStore=time.tStore; % when to save
        satStore=zeros(length(solver.nodes),length(tStore));
        headStore=zeros(size(satStore));
    else
        % when to save
        if length(time.tStore)==1
            % store every tStore until the end
            tStore=time.tStore:time.tStore:time.tStop;
        else
            % tStore is a vector of stores
            tStore=time.tStore;
        end
        headStore=zeros(geom.nz,length(tStore));
        satStore=zeros(geom.nz,length(tStore));
        timestore=zeros(1,length(tStore));
    end    
end

% Get a vector of times that have to be hitted exactly
mustDos=sort([boundary.changeTimes(2:end),tStore,time.tStop]);
a=mustDos(2:end)-mustDos(1:end-1);
a= a~=0;
mustDos=mustDos(a);
if isempty(mustDos)
    mustDos=time.tStop;
end

% Calcuate the hamonic average of the saturated permeabilities
Ksat=zeros(geom.nz+1,1);
Ksat(2:geom.nz)=2*param.kin(2:geom.nz).*param.kin(1:geom.nz-1)./(param.kin(2:geom.nz)+param.kin(1:geom.nz-1));
% bottom and top are the originals
Ksat(1)=param.kin(1);Ksat(end)=param.kin(end);

% Initialize boundaries
boundary.BottomHead=boundary.headBottom(1);
boundary.TopFlux=boundary.fluxTop(1);
boundary.kr=0;

% and the timestep
dt=time.dt(1);
% and check which to update during the run
if length(boundary.headBottom)>1;   doBottom=true;
else                                doBottom=false; end
if length(boundary.fluxTop)>1;      doTop=true;
else                                doTop=false;    end
if length(time.dt)>1;               doTime=true;
else                                doTime=false;   end


%% --------------------------Time loop --------------------------------
while tim<time.tStop
    
        
    %% initialze, print and make the system ready
    
    % Check what boundry condtion to use at the bottom
    if any(boundary.changeTimes==tim)  || tscount==1
        if doBottom
            tmp=boundary.headBottom(boundary.changeTimes<=tim);
            boundary.BottomHead=tmp(end);
        end
        if doTop
            tmp=boundary.fluxTop(boundary.changeTimes<=tim);
            boundary.TopFlux=tmp(end);
            
            if boundary.doAdHoc1LayerRoots
                % TEMPORARY; FINALIZE OR REMOVE
                tmp=boundary.adHocUptakeSerie(boundary.changeTimes<=tim);
               boundary.adHocUptake=tmp(end);                
            end
            
        end
        if doTime
            tmp=time.dt(boundary.changeTimes<=tim);
            dt=tmp(end);
        end
        
        % as the boundary Krel does not change, calculate it once here
        if ~boundary.doFreeDrain && strcmp(param.what,'VG')
            aB=param.alpha(1).*abs(boundary.BottomHead);
            bB=(1+aB^param.n(1)).^param.m(1);
            boundary.kr=((1 -((aB)^(param.n(1)-1))/bB)^2)*(bB^param.aKR(1));
        elseif ~boundary.doFreeDrain && strcmp(param.what,'RG')
            boundary.kr=exp(-param.alpha(1).*abs(boundary.BottomHead));          
        end
        
        
        % When the boundaries are changed, recalculate the old rel. perm.
        old.kr = newtonKrel(headinter,se,param,boundary);
    end
    
    
    
    % Check for the maximum possible time step size
    a=mustDos-tim;
    a(a<=0)=inf;
    dtMax=min(a);
    dt(1.05*dt>dtMax)=dtMax;
    
    if doPrint
        % print out the info ever X tscount
        if mod(tscount,100)==0; fprintf(1,[num2str(round(100*100*tim/time.tStop)/100) ' '])
            if prcount == 5; fprintf(1,'\n'); prcount=0;
            else prcount=prcount+1; end
        end
    end
    
    % Set the water content at the old timestep
    wc_old=wc_new;
    old.wc=wc_old;
    
    %% ----------------------- Start of Nonlinear Iteration --------
    
    % Initialize a counter
    zaehler1=0;
    zaehler2=0;
    
    %-------------------------------------------------------
    % While loop, do as many Newton-Raphson steps until convercence is reached.
    % Check every solver.convBreak1 steps if convergance problems and break
    % out if solver.convBreak2 is reached
    % The Newton-Raphson Iterations are done to find the minimum for
    % (S Storage h^n+1-h^n)/dt + nf S^n+1-S^n/dt) - div flux^n+1) =
    % f(h^n+1).
    
    % Initialize some variables
    rhs=1;
    rhsFirst=1;
    incrOrg=1;
    old.rhs=[];
    doBreak=false;
    
    
    % Keep iterating until the defect between current and last iteration
    % step is small enough. Second break of criteria is present only aviod
    % a crash when reaching steady state (hence no defect decrease)
    while sqrt(sum(rhs.^2))/sqrt(sum(rhsFirst.^2))>solver.eps && max(abs(incrOrg))>solver.eps/10 && sqrt(sum(rhs.^2)) > solver.eps/10
        %max(abs(incrOrg))>solver.eps/100%
        
        %%  Do the iteration steps
        
        % Calculate the Jacobian and the right hand side for the Newton iterations
        [J,rhs]=newtonJacobian(headinter,Ksat,wc_old,delzNN,delzCS,dt,boundary,param,solver.epsh,head_old,old);
        
        % solve for the increment
        incr=J\(-rhs);
        incrOrg=incr;
        
        [~, msgid] = lastwarn;
        lastwarn('')
        
        if any(isinf(incr)) || any(isnan(incr)) || strcmp(msgid,'MATLAB:singularMatrix') || strcmp(msgid,'MATLAB:nearlySingularMatrix')
             % Jacobian is singular
            if  sqrt(sum(rhs.^2)) < solver.eps/10
                % clean anyway
                break
            end
             
             if solver.jacSingBreak
               % Break out!
                doBreak=true;
                if solver.doEnKF
                    headinter=ones(size(headinter))*inf;
                elseif solver.doDream
                    satStore=ones(size(satStore))*inf;                    
                else
                    headStore=ones(size(headStore))*inf;
                    satStore=ones(size(satStore))*inf;
                end        
                
                disp('---------------Jacobian singular, breaking out!')
                break
                
            else               
                % restart (i.e. make the system go to restart)
                incr=zeros(geom.nz,1);
                rhs=ones(geom.nz,1)*1e-20;
                zaehler1=solver.convBreak1-1;
            end
        end
        
        % set the first defect to the current if first iteration step
        if zaehler1 == 0; rhsFirst=rhs; end
        
      
           
        
        %% Line search
        % decrease the increment if the solution deteriorates
        % implemeted after the numerics script of O. Ippisch p.140
        if solver.doLineSearch           
            % initialize alpha
            alpha=1;
            
            % original norm
            normf0=norm(rhs,inf);
            
            % alpha=1 norm
            [sat,se]=calcSat((headinter+incr),param);
            wc_new=sat.*param.poro;
            storeX=store((headinter+incr),wc_new,wc_old,param,dt,delzCS,head_old);
            kr = newtonKrel((headinter+incr),se,param,boundary);
            rhs=richardsNewton((headinter+incr),kr,Ksat,storeX,delzNN,boundary);            
            normf1=norm(rhs,inf);
            
            % line search loop
            while normf1 >= normf0 && alpha >= 2*1/10000 && normf1>1e-5
                alpha=alpha*0.5;
                [sat,se]=calcSat((headinter+alpha*incr),param);
                wc_new=sat.*param.poro;
                storeX=store((headinter+alpha*incr),wc_new,wc_old,param,dt,delzCS,head_old);
                kr = newtonKrel((headinter+alpha*incr),se,param,boundary);
                rhs=richardsNewton((headinter+alpha*incr),kr,Ksat,storeX,delzNN,boundary);               
                normf1=norm(rhs,inf);
            end
            
            % decrease the increment
            incr=alpha*incr;
        else
            % if no line search, just calculate the defect
            [sat,se]=calcSat((headinter+incr),param);
            wc_new=sat.*param.poro;
            storeX=store((headinter+incr),wc_new,wc_old,param,dt,delzCS,head_old);
            kr = newtonKrel((headinter+incr),se,param,boundary);
            rhs=richardsNewton((headinter+incr),kr,Ksat,storeX,delzNN,boundary);
        end
        
        
        %% back to the original solver        
              
        % Save same variables for later re-use
        old.rhs=rhs;
        old.kr=kr;
        old.store=storeX;
        
        % Update the head
      %  hx=headinter;
        headinter=headinter+incr;
    %    incrOLD=incrOrg;
        
        % Increase the counter by one
        zaehler1=zaehler1+1;
        
        % Check all X iteration for convergence problem. If such, decrease
        % timestep size and restart the timestep  
        
        if(mod(zaehler1,solver.convBreak1)==0)
           
            % reduce with 60%
            if dt~=time.dtMin
                dt=max(time.dtMin,0.4*dt);      
                % reset the counters
                zaehler1=0;
            end
            
            zaehler2=zaehler2+zaehler1;
            if zaehler2 >= solver.convBreak2
                % Quite and bail out               
                if solver.doEnKF
                    headinter=ones(size(headinter))*inf;                
                else
                    headStore=ones(size(headStore))*inf;
                    satStore=ones(size(satStore))*inf;
                end        
                disp('---------------Conv. probl., breaking out!')
                
                doBreak=true;
                break
            end
            
            
            % otherwise restart and restart
            headinter=head_old;
            [sat,se]=calcSat(headinter,param);
            clear rhs kr storeX
            rhs=1000000;
            wc_new=sat.*param.poro;
            wc_old=wc_new;
            old.wc=wc_new;
            old.store=store(headinter,wc_new,wc_old,param,dt,delzCS,head_old);
            old.kr = newtonKrel(headinter,se,param,boundary);
            old.rhs=[];
            
            if doPrint
                disp(['---------------Conv.probl., new dt: ' num2str(dt)])
            end
        end
        
        %-------------- end of Newton iteration ---------------
    end
    %%
    
    if doBreak
        break
    end
    
    % Set time to the new time
    tim=tim+dt;
    
    zTot=zTot+zaehler1+zaehler2;
    if zTot>time.MaxIter && tim<time.tStop
        % break off!
        if solver.doEnKF
            headinter=ones(size(headinter))*inf;
        elseif solver.doDream
            satStore=ones(size(satStore))*inf;
            headStore=ones(size(headStore))*inf;
        else
            headStore=ones(size(headStore))*inf;
            satStore=ones(size(satStore))*inf;
        end
        disp('---------------Too many iterations, breaking out')
        break
    end
    
%     if feldcount>=190 && tim >= tStore(feldcount)
%         keyboard
%     end
    
    % FOR COUPLING: CALCULATE THE FLUX FROM THE SECOND LAST TO THE LAST
    % CELL AND USE THIS AS A REPLACEMENT FOR A RECHARGE
      % flux
%       fx=kr(2)*Ksat(1)*((headinter(2)-headinter(1))./delzNN(2)+1);
%     R=R+fx*dt;


% version 2 (new edition 20.02.2018 to better account for backwards flow
if headinter(1)-boundary.headBottom>-delzNN(1); krR=kr(2); else, krR=kr(1); end
     fx=krR*Ksat(1)*((headinter(1)-boundary.headBottom)./delzNN(1)+1);
     R=R+fx*dt;


% version 3: better mass balance on shallow columns with a ghost cell?? (01.03.2017)
%      fx=kr(1)*Ksat(1)*((headinter(1)-(boundary.headBottom+delzNN(1)))./delzCS(1)+1);
%         R=R+fx*dt;

        
    %  -------------Output after given time -------------------------
    if ~solver.doEnKF && tim >= tStore(feldcount)
        
        % Saturation renamed
        if solver.doDream
            % Saturation renamed
            satStore(:,feldcount)=sat(solver.nodes);
            headStore(:,feldcount)=headinter(solver.nodes);
        else
            
            % Head renamed
            headStore(:,feldcount)=headinter;
            
            % Saturation renamed
            satStore(:,feldcount)=sat;
            
            % timestore
            timestore(feldcount)=tim;
        end
        
        % Set outputcounter one up
        feldcount=feldcount+1;
    end
    
    % Set time counter one up
    tscount = tscount + 1;
    
    %Check if the next timestep should be different
    if time.doAdaptiveTime
        dt=timecheck(dt,zaehler1,sat-wc_old./param.poro);
        % and check the boundaries
        dt=min(time.dtMax,dt);
        dt=max(time.dtMin,dt);
    end
    
    % Set old time to the actual time to start the new time.
    head_old=headinter;
    
    
    %----------------------- End of time loop -----------------------------
end


if solver.doEnKF
    % hence if there is only suppose to be one output at the end of the
    % time
    headStore=headinter;
end