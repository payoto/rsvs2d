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
        =ExecuteSnakes_Optim(entryPoint,gridrefined,looprestart,baseGrid,connectstructinfo...
        ,param,paramspline,outinfo,nIter,nProf,nPop)
    
    
    switch entryPoint
        case 'snak'
            [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
                =ExecuteSnakes_Optim_snak(gridrefined,looprestart,baseGrid,connectstructinfo...
                ,param,paramspline,outinfo,nIter,nProf,nPop);
        case 'sens'
            [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
                =ExecuteSnakes_Optim_sens(gridrefined,looprestart,baseGrid,connectstructinfo...
                ,param,paramspline,outinfo,nIter,nProf,nPop);
            
    end
    
    
end



function [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
        =ExecuteSnakes_Optim_snak(gridrefined,looprestart,baseGrid,connectstructinfo...
        ,param,paramspline,outinfo,nIter,nProf,nPop)
    % Executes the snakes edge detection process
    %
    
    procStr='SNAKE PROCESS';
    
    [textOut1,tStart]=evalc('PrintStart(procStr,2);');
    
    varExtract={'refineSteps','snakData'};
    [refineSteps,snakData]=ExtractVariables(varExtract,param);
    
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
    tecStruct.snakposition=snakposition;
    tecStruct.baseGrid=baseGrid;
    tecStruct.nPop=nPop;
    tecStruct.fineGrid=gridrefined;
    
    if strcmp('all',snakData)
        [outinfo,snakSave]=FullResultsRequest(gridrefined,connectstructinfo,baseGrid,...
            snakSave,param,outinfo,nIter,nProf,looprestart,restartsnake,tecStruct);
    else
        outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
            restartsnake,snakSave,tecStruct);
    end
    [textOut2,~]=evalc('PrintEnd(procStr,2,tStart)');
    fprintf([textOut1,textOut,textOut2])
end




function [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
        =ExecuteSnakes_Optim_sens(gridrefined,looprestart,baseGrid,connectstructinfo...
        ,param,paramspline,outinfo,nIter,nProf,nPop)
    % Executes the snakes edge detection process
    %
    
    procStr='SNAKE PROCESS';
    
    [textOut1,tStart]=evalc('PrintStart(procStr,2);');
    
    varExtract={'refineSteps','snakData'};
    [refineSteps,snakData]=ExtractVariables(varExtract,param);
    
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
    tecStruct.snakposition=snakposition;
    tecStruct.baseGrid=baseGrid;
    tecStruct.nPop=nPop;
    tecStruct.fineGrid=gridrefined;
    
    if strcmp('all',snakData)
        warning('Incompatible parameters have been declared: snakData request ''all'' ignored')
        outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
            restartsnake,snakSave,tecStruct);
    else
        outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
            restartsnake,snakSave,tecStruct);
    end
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
    varExtract={'typeBound','subdivType','resampleSnak','TEShrink','LEShrink','typeCorner'};
    [typeBound,subdivType,resampleSnak,TEShrink,LEShrink,typeCorner]=ExtractVariables(varExtract,param);
    if ~exist('typeBound','var'), typeBound='vertex'; end
    sharpen=[LEShrink,TEShrink];
    
    for ii=1:length(loop)
        startPoints=loop(ii).(typeBound).coord;
        loop(ii).isccw=CCWLoop(startPoints);
        [newPoints,projPoints]=SubDivision(startPoints,refineSteps,subdivType,sharpen,typeCorner);
        loop(ii).subdivision=newPoints;
        loop(ii).subdivspline=newPoints;
        if resampleSnak
            resampPoints=ResampleSpline(projPoints,paramspline);
            loop(ii).subdivspline=resampPoints;
        end
        %newPoints=SubSurfBSpline(startPoints(1:end-2,:),refineSteps);
        
    end
end


%% Handle full results request

function [outinfo,snakSave]=FullResultsRequest(gridrefined,connectstructinfo,baseGrid,...
        snakSave,param,outinfo,nIter,nProf,looprestart,restartsnake,tecStruct)
    snakSave2=snakSave;
    volfra=snakSave(end).volumefraction;
    snakSave(end).volumefraction=struct('targetfill',[],'currentfraction'...
        ,[],'totVolume',[]);
    snakSave(end).volumefraction.targetfill=[volfra(:).targetfill];
    snakSave(end).volumefraction.currentfraction=[volfra(:).volumefraction];
    snakSave(end).volumefraction.totVolume=[volfra(:).totalvolume];
    
    outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
        restartsnake,snakSave,tecStruct);
    
    
    writeDirectory=outinfo.dirprofile;
    marker=[int2str(nIter),'_',int2str(nProf)];
    
    % TecPlot Data
    [fidTecPLT,pltFileName]=OpenTecPLTFile(writeDirectory,marker);
    fidTecLAY=OpenTecLayFile(writeDirectory,marker);
    PersnaliseLayFile(fidTecLAY,pltFileName)
    
    TecplotOutput('snakes',fidTecPLT,baseGrid,gridrefined,snakSave2,connectstructinfo)
    
    
end
