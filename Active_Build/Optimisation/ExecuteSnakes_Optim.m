%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%       Optimisation Using 
%          Parametric Snakes for
%      for Aerodynamic shape
%         parametrisation
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
        =ExecuteSnakes_Optim(gridrefined,looprestart,baseGrid,connectstructinfo,param,paramspline,outinfo,nIter,nProf)
    % Executes the snakes edge detection process
    %
    
    procStr='SNAKE PROCESS';
    
    [textOut1,tStart]=evalc('PrintStart(procStr,2);');
    
    varExtract={'refineSteps'};
    [refineSteps]=ExtractVariables(varExtract,param);
    
    callerString='Snakes(gridrefined,looprestart,baseGrid,connectstructinfo,param);';
    [textOut,snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=evalc(callerString);
    
    if length(loopsnaxel)==length(looprestart)
        for ii=1:length(loopsnaxel)
            looprestart(ii).snaxel=loopsnaxel(ii).snaxel;
        end
    else
        looprestart=loopsnaxel;
    end
    
    
    looprestart=SubdivisionSurface_Snakes(looprestart,refineSteps,param,paramspline);
    
    outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,restartsnake,snakSave);
    [textOut2,~]=evalc('PrintEnd(procStr,2,tStart)');
    fprintf([textOut1,textOut,textOut2])
end


%% Print to screen functions

function [tStart]=PrintStart(procStr,lvl)
    
    procStart=[procStr,' start'];
    tStart=now;
    switch lvl
        case 1
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp('--------------------------------------------------------------------------------------------')
            disp(procStart)
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
            
        case 2
            disp('  ')
            disp('-----------------------')
            disp(procStart)
            
        case 3
            disp('----------')
            disp(procStart)
    end
    
end

function [tElapsed]=PrintEnd(procStr,lvl,tStart)
    
    procStart=[procStr,' end'];
    tEnd=now;
    tElapsed=tEnd-tStart;
    switch lvl
        case 1
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(procStart)
            disp('--------------------------------------------------------------------------------------------')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
            
        case 2
            
            
            disp(['Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(procStart)
            disp('-----------------------')
            disp('  ')
        case 3
            disp('----------')
            
    end
    
end

%% Subdivision process

function [loop]=SubdivisionSurface_Snakes(loop,refineSteps,param,paramspline)
    % Function taking in a closed loop of vertices and applying the subdivision
    % process
    % typBOund is te type of boundary that is expected, it can either be the
    % string 'vertex' (default) or the string 'snaxel' to show that the snaxel
    % process has been used
     varExtract={'typeBound','subdivType'};
    [typeBound,subdivType]=ExtractVariables(varExtract,param);
    if ~exist('typeBound','var'), typeBound='vertex'; end
    
    for ii=1:length(loop)
        startPoints=loop(ii).(typeBound).coord;
        loop(ii).isccw=CCWLoop(startPoints);
        [newPoints,projPoints]=SubDivision(startPoints,refineSteps,subdivType);
        resampPoints=ResampleSpline(projPoints,paramspline);
        %newPoints=SubSurfBSpline(startPoints(1:end-2,:),refineSteps);
        loop(ii).subdivision=newPoints;
        loop(ii).subdivspline=resampPoints;
    end
end
