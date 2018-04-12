% Part of a new version of the solution for the Richards equation
% Written Dec 2011 by Daniel Erdal

% Version 1.0
% Last Documented edition: 30.05.12 (DE)


function dt=timecheck(dtOld,counter,satDiff)

%%
% With current setup aimed at the newton iterations
% (with picard, 100 iterations is not so much...)


if counter > 100 % || max(abs(satDiff)) > 0.1
    % lower the time step size
    if counter > 2000
         dt=0.62*dtOld;
    elseif counter > 1000
         dt=0.74*dtOld;
    elseif counter > 500
         dt=0.86*dtOld;
    elseif counter > 250
         dt=0.98*dtOld;
    else
         dt=0.99*dtOld;
    end      
   
elseif counter < 50
    % increase the time step size
    if counter >40
        dt=1.001*dtOld;
    elseif counter > 30
        dt=1.002*dtOld;
    elseif counter > 20
        dt=1.003*dtOld;
    elseif counter > 10
        dt=1.005*dtOld;
    elseif counter > 4
        dt=1.007*dtOld;
    else
       dt=1.01*dtOld;
    end
else
    dt=dtOld;
end



