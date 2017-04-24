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
        case 'post'
            [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
                =ExecuteSnakes_Optim_post(gridrefined,looprestart,baseGrid,connectstructinfo...
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
    tecStruct.volumefraction=snakSave(end).volumefraction;
    tecStruct.snaxel=restartsnake.snaxel;
    
    if strcmp('all',snakData)
        [outinfo,snakSave]=FullResultsRequest(gridrefined,connectstructinfo,baseGrid,...
            snakSave,param,outinfo,nIter,nProf,looprestart,restartsnake,tecStruct);
    else
        
        outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
            restartsnake,snakSave,tecStruct);
    end
    [textOut2,~]=evalc('PrintEnd(procStr,2,tStart)');
    disp([textOut1,textOut,textOut2])
end


function [tecsnaxel,tecsnakposition,snakSave,looprestart,restartsnake,outinfo]...
        =ExecuteSnakes_Optim_sens(gridrefined,supportstructsens,baseGrid,connectstructinfo...
        ,param,paramspline,outinfo,nIter,nProf,nPop)
    % Executes the snakes edge detection process
    %
    % loop
    
    procStr='SNAKE PROCESS';
    
    [textOut1,tStart]=evalc('PrintStart(procStr,2);');
    
    varExtract={'refineSteps','snakData','axisRatio'};
    [refineSteps,snakData,axisRatio]=ExtractVariables(varExtract,param);
    
    looprestart=supportstructsens.loopsens;
    restartsnake=supportstructsens;
    
    % callerString='Snakes(gridrefined,looprestart,baseGrid,connectstructinfo,param);';
    % [textOut,snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=evalc(callerString);
    
    
    [looprestart]=FinishLoops(looprestart,param);
    looprestart=SubdivisionSurface_Snakes(looprestart,refineSteps,param,paramspline);
    
    
    [tecsnaxel,tecsnakposition]=LooptoTecSnax(looprestart);
    
    tecStruct.snakposition=tecsnakposition;
    tecStruct.baseGrid=baseGrid;
    tecStruct.nPop=nPop;
    tecStruct.fineGrid=gridrefined;
    tecStruct.volumefraction.targetfill=...
        [supportstructsens.volumefraction(:).targetfill];
    tecStruct.volumefraction.currentfraction=...
        [supportstructsens.volumefraction(:).volumefraction];
    tecStruct.volumefraction.totVolume=...
        [supportstructsens.volumefraction(:).totalvolume];
    tecStruct.snaxel=tecsnaxel;
    
    snakSave=struct('tecstruct',tecStruct,'warn','This is not a standard snakSave');
    
    if strcmp('all',snakData)
        warning('Incompatible parameters have been declared: snakData request ''all'' ignored')
        outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
            restartsnake,snakSave,tecStruct);
    else
        outinfo=OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,...
            restartsnake,snakSave,tecStruct);
    end
    [textOut2,~]=evalc('PrintEnd(procStr,2,tStart)');
    fprintf([textOut1,'\n   Sensitivity Profile post-treated\n',textOut2])
end

% function [looprestart]=ApplyAxisRatioToSens(looprestart,axisRatio)
%     
%      for ii=1:numel(loop)
%          
%          
%          
%      end
% end

function [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
        =ExecuteSnakes_Optim_post(gridrefined,looprestart,baseGrid,connectstructinfo...
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
    tecStruct.volumefraction=snakSave(end).volumefraction;
    tecStruct.snaxel=restartsnake.snaxel;
    
    if strcmp('all',snakData)
        [outinfo,snakSave]=FullResultsRequest(gridrefined,connectstructinfo,baseGrid,...
            snakSave,param,outinfo,nIter,nProf,looprestart,restartsnake,tecStruct);
    else
        
        outinfo=OptimisationOutput('profile',param,outinfo,looprestart,...
            restartsnake,snakSave,tecStruct);
    end
    [textOut2,~]=evalc('PrintEnd(procStr,2,tStart)');
    disp([textOut1,textOut,textOut2])
end


function [tecsnaxel,tecsnakposition]=LooptoTecSnax(loopsens)
    tot=0;
    for ii=1:length(loopsens)
        tot=tot+numel(loopsens(ii).snaxel.index);
    end
    
    tecsnaxel=repmat(struct('index',[],'snaxnext',[],'v',1e-10),[1 tot]);
    tecsnakposition=repmat(struct('index',[],'coord',[],'vector',[]),[1 tot]);
    
    kk=1;
    for ii=1:length(loopsens)
        nLoop=numel(loopsens(ii).snaxel.index);
        for jj=1:numel(loopsens(ii).snaxel.index)
            
            tecsnaxel(kk).index=loopsens(ii).snaxel.index(jj);
%             nextSub=[kk+1:min([kk+1,nLoop]),mod(kk,nLoop)+1:1];
            tecsnaxel(kk).snaxnext=loopsens(ii).snaxel.snaxnext(jj);
            %tecsnaxel(kk).snaxnext=loopsens(ii).snaxel.index(nextSub);
            
            tecsnakposition(kk).index=loopsens(ii).snaxel.index(jj);
            tecsnakposition(kk).coord=loopsens(ii).snaxel.coordnoscale(jj,:);
            tecsnakposition(kk).vector=loopsens(ii).snaxel.vector(jj,:);
            kk=kk+1;
        end
    end
    
    
    
    
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
        %startPoints=loop(ii).(typeBound).coord;
        isccw=CCWLoop(loop(ii).(typeBound).coord);
        if ~isccw
           loop(ii).(typeBound).coord=flip(loop(ii).(typeBound).coord);
           isccw=CCWLoop(loop(ii).(typeBound).coord);
        end
        loop(ii).isccw=isccw;
        startPoints=loop(ii).(typeBound).coord;
        [newPoints,projPoints]=SubDivision(startPoints,refineSteps,subdivType,sharpen,typeCorner);
        if isempty(newPoints)
            disp('Loop is a single point')
            newPoints=startPoints;
            projPoints=startPoints;
        end
        loop(ii).subdivision=newPoints;
        loop(ii).subdivspline=projPoints;
        if resampleSnak
            resampPoints=ResampleSpline(projPoints,paramspline);
            % DVP Check for resampling if using mesh motion.
            loop(ii).subdivspline=resampPoints;
        end
        %newPoints=SubSurfBSpline(startPoints(1:end-2,:),refineSteps);
        
    end
    if ~all(~sharpen)
        %
        rmLoop=false(size(loop));
        for ii=1:numel(loop)
            rmLoop(ii)=isempty(loop(ii).isccw);
        end
        loop(rmLoop)=[];
        if any(rmLoop)
            warning(' %i empty loops removed.',sum(rmLoop))
        end
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
    
    outinfo=OptimisationOutput('profilepost',param,outinfo,looprestart,...
        restartsnake,snakSave,tecStruct);
    
    
    writeDirectory=outinfo.dirprofile;
    marker=[int2str(nIter),'_',int2str(nProf)];
    
    % TecPlot Data
    [fidTecPLT,pltFileName]=OpenTecPLTFile(writeDirectory,marker);
    fidTecLAY=OpenTecLayFile(writeDirectory,marker);
    PersnaliseLayFile(fidTecLAY,pltFileName)
    
    TecplotOutput('snakes',fidTecPLT,baseGrid,gridrefined,snakSave2,connectstructinfo)
    
    
end
