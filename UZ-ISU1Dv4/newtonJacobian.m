% Part of a new version of the solution for the Richards equation
% Written Dec 2011 by Daniel Erdal
% Compacted version written May 2012

% Version 2.0
% Last Documented edition: 05.12.13 (DE)

% Calculate the entries for the Jacobian
% assuming knoledge of entries on the main and 1st off diagonals

% Output is the Jacobian and the right hand side vector

function [J,rhs]=newtonJacobian(h,Ksat,wc_old,delzNN,delzCS,dt,boundary,param,epsh,h_old,old)
%% 
 global doFlow
    
if boundary.doFreeDrain && boundary.freeDrainOption == 3
    
    % To assure that the bottom does not have different drainge in the
    % calculation of the jacobian is it required that the 0.9 saturation is
    % fullfilled from the start
    if  strcmp(param.what,'VG')    
        s1=1./((1+(boundary.alphaB*abs(h(1))).^boundary.nB).^boundary.mB);
    elseif  strcmp(param.what,'RG')
         s1=(exp(-0.5*boundary.alphaB*abs(h(1))).*(1+0.5*boundary.alphaB*abs(h(1)))).^(2/2.5);
    end
    s1(s1>1)=1;
    s1(h(1)>0)=1;
    
    if s1>0.90
       doFlow=true;
    else
        doFlow=false;
    end
end


nz=length(h);

% During iterations, reuse the values calculated when evaluating the defect
% in the "solverStatic" script. If within iterations, only reuse water
% content and permeabilities.
if isempty(old.rhs)
    wc_new=old.wc;
    kr_new=old.kr;
    store_new=store(h,wc_new,wc_old,param,dt,delzCS,h_old); 
    f0=richardsNewton(h,kr_new,Ksat,store_new,delzNN,boundary);
else  
    kr_new=old.kr;
    store_new=old.store;
    f0=old.rhs;
end

% setup indecies
a=1:3:nz;b=2:3:nz;c=3:3:nz;

% add the delta h
delta=-epsh*h-epsh+5*eps;
delta(abs(delta)<100*eps) = epsh;
hx=h+delta;
ha=h; ha(a)=hx(a);
hb=h; hb(b)=hx(b);
hc=h; hc(c)=hx(c);




% Calculate the new saturations and water contents for the head change
[sat,se]=calcSat(hx,param);
wc_new_x=param.poro.*sat;

% Storage and relative permeability for the increased head
store_new_x=store(hx,wc_new_x,wc_old,param,dt,delzCS,h_old);
kr_new_x = newtonKrel(hx,se,param,boundary);

% Combine new and old values to create three cases (3 cases needed due to
% offdiagonals that otherwise becomes wrong)
kr_new_a=kr_new; kr_new_a(a+1)=kr_new_x(a+1); 
kr_new_b=kr_new; kr_new_b(b+1)=kr_new_x(b+1); 
kr_new_c=kr_new; kr_new_c(c+1)=kr_new_x(c+1); 

store_new_a=store_new; store_new_a(a)=store_new_x(a); 
store_new_b=store_new; store_new_b(b)=store_new_x(b); 
store_new_c=store_new; store_new_c(c)=store_new_x(c); 

% Caluculate the function f for three cases (3 cases needed due to
% offdiagonals that otherwise becomes wrong)
xa=richardsNewton(ha,kr_new_a,Ksat,store_new_a,delzNN,boundary);
xb=richardsNewton(hb,kr_new_b,Ksat,store_new_b,delzNN,boundary);
xc=richardsNewton(hc,kr_new_c,Ksat,store_new_c,delzNN,boundary);

%% Calculate the derivatives
% diagonal (fi(h+delta(i))-fi(h))/delta(i)
fp=zeros(nz,1);
fp(a)=xa(a);
fp(b)=xb(b);
fp(c)=xc(c);

dia=(fp-f0)./delta;

% offdiagonal +1; (fi(h+delta(i+1))-fi(h))/delta(j)
fp=zeros(nz-1,1);
fp(1:3:nz-1)=xb(1:3:nz-1);
fp(2:3:nz-1)=xc(2:3:nz-1);
fp(3:3:nz-1)=xa(3:3:nz-1);

diaOff1=(fp-f0(1:end-1))./delta(2:end);

% offdiagonal -1; (fj(h+delta(i))-fj(h))/delta(i)
fp=zeros(nz-1,1);
fp(1:3:nz-1)=xa(2:3:nz);
fp(2:3:nz-1)=xb(3:3:nz);
fp(3:3:nz-1)=xc(4:3:nz);

diaOff2=(fp-f0(2:end))./delta(1:end-1);

%% Assemble the Jacobian
J=sparse([1:nz,2:nz,1:nz-1],[1:nz,1:nz-1,2:nz],[dia;diaOff2;diaOff1]);
rhs=f0;


