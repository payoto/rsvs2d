 

function [varargout]=SnakCutCell_SurrogateHandler(entryPoint,varargin)
    
    P = mfilename('fullpath');
    lengthName=length('SnakCutCell_SurrogateHandler');
    P(end-lengthName:end)='';
    P=[P,'\..\..\'];
    startDir=cd(P);
    
    switch entryPoint
        case 'init'
            % paramsur,outinfo
            [varargout{1}]=SnakCutCell_init(varargin{:});
        case 'run'
            % datastruct,parsurrogate,supportstructure,outinfo
            [varargout{1}]=SnakCutCell_Run(varargin{:});
        otherwise
            error('Undefined output case')
    end
    
    cd(startDir)
end

%% Initialisation

function [supportstructure]=SnakCutCell_init(paramsur,outinfo)
    
    % Initialise Workspace
    include_EdgeInformation
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    include_PostProcessing
    include_Mex_Wrapper
    include_Optimisation
    
    varExtract={'snakCase','worker'};
    [snakCase,worker]=ExtractVariables(varExtract,paramsur);
    paramsnak=structInputVar(snakCase);
    paramsnak.initparam=DefaultSnakeInit(paramsnak);
    paramsnak.spline=DefaultOptimSpline();
    % Start Parrallel pool
    if numel(gcp('nocreate'))==0
        comStr=computer;
        if strcmp(comStr(1:2),'PC')
            poolName=parallel.importProfile('ExportOptimSnakes.settings');
        else
            poolName=parallel.importProfile('ExportOptimSnakesLinux.settings');
        end
        clusterObj=parcluster(poolName);
        clusterObj.NumWorkers=worker;
        saveProfile(clusterObj);
        parpool(poolName)
    end
    
    % Initialisation structure
    [unstructured,baseGrid,gridrefined,connectstructinfo,...
        unstructuredrefined,loop]=GridInitRefine(paramsnak);
    [~,~,~,~,restartsnake]=ExecuteSnakes_Surrogate(gridrefined,loop,...
        baseGrid,connectstructinfo,paramsnak.initparam,paramsnak.spline,outinfo,0);
    
    % Support Structure
    supportstructure.gridrefined=gridrefined;
    supportstructure.restartsnake=restartsnake;
    supportstructure.baseGrid=baseGrid;
    supportstructure.connectstructinfo=connectstructinfo;
    supportstructure.param=paramsnak;
end

function [unstructured,unstructReshape,gridrefined,connectstructinfo,...
        unstructuredrefined,loop]=GridInitRefine(param)
    % Grid initialisation
    
    [unstructured,~,unstructReshape]=...
        GridInitialisationV2(param);
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    
    
    
end

function [paraminit]=DefaultSnakeInit(paraminit)
    
    paraminit.snakes.step.snakesSteps=40;
    paraminit.general.restart=false;
end

function [paroptimspline]=DefaultOptimSpline()
    
    
    paroptimspline.splineCase='aerosnake';
    paroptimspline.domain='normalizeX';
end
%% Running

function [datastruct]=SnakCutCell_Run(datastruct,parsurrogate,supportstructure,outinfo,varargin)
    
    indList=find([datastruct(:).flag]<0);
    [datastruct(indList)]=PerformIteration(parsurrogate,outinfo,...
        datastruct(indList),indList,supportstructure);
    indList=indList(([datastruct(indList).flag]<=1));
    
    if numel(varargin)==1
        [datastruct(indList).flag]=deal(varargin{1});
    end
    
end

function [runcases]=PerformIteration(parsurrogate,outinfo,runcases,indList,supportstructure) %gridrefined,restartsnake,...
        %baseGrid,connectstructinfo)
    
    
    nPop=length(runcases);
    
    paramsnake=supportstructure.param;
    paramspline=paramsnake.spline;
    gridrefined=supportstructure.gridrefined;
    restartsnake=supportstructure.restartsnake;
    baseGrid=supportstructure.baseGrid;
    connectstructinfo=supportstructure.connectstructinfo;
    [captureErrors{1:nPop}]=deal('');
    
    parfor ii=1:nPop
        %for ii=1:nPop
        
        currentMember=[runcases(ii).fill,0.5 0.35 0.2 0.05 1e-4];
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        try
            % Normal Execution
            runcases(ii)=NormalExecutionIteration(runcases(ii),newRefGrid,newrestartsnake,...
                newGrid,connectstructinfo,paramsnake,paramspline,outinfo,indList(ii),parsurrogate);
            runcases(ii).flag=1;
            
        catch MEexception
            % Error Capture
            runcases(ii).constraint=0;
            runcases(ii).exception=['error: ',MEexception.identifier];
            captureErrors{ii}=MEexception.getReport;
        end
    end
    % Add Invalidity
    dataMat=[runcases(:).additional];
    resPos=dataMat(:,2)>0;
    [runcases(resPos).flag]=deal(0);
    
    
end

function [newGrid,newRefGrid,newRestart]=ReFillGrids(baseGrid,refinedGrid,restartsnake,connectstructinfo,newFill)
    
    activeCell=logical([baseGrid.cell(:).isactive]);
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    cellCentreInd=[restartsnake.cellCentredGrid(:).index];
    
    if numel(newFill)~=numel(activeCellSub)
        error('Fill and Active Set do not match in size')
    end
    if sum(abs(newFill))==0
        newFill(round(numel(newFill)/2))=1e-3;
    end
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    newRestart=restartsnake;
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFill(ii);
        newRestart.volfracconnec.cell(activeCellSub(ii)).targetfill=newFill(ii);
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFill(ii));
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],cellCentreInd);
        [newRestart.cellCentredGrid(newSub).fill]=deal(newFill(ii));
    end
    
end

function population=NormalExecutionIteration(population,newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,ii,paramsurrogate)
    
    varExtract={'restart','boundstr'};
    [isRestart,boundstr]=ExtractVariables(varExtract,paramsnake);
    
    if ~isRestart
        [newRefGrid]=EdgePropertiesReshape(newRefGrid);
        [newGrid]=EdgePropertiesReshape(newGrid);
        [newrestartsnake]=GenerateEdgeLoop(newRefGrid,boundstr,true);
    end
    
    [snaxel,snakposition,snakSave,loop,restartsnake,outTemp]=...
        ExecuteSnakes_Surrogate(newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,ii);
    population.location=outTemp.rootDir;
    
    
    [objValue,additional]=CutCellFlow(paramsurrogate,population,loop);
    additional.snaxelVolRes=snakSave(end).currentConvVolume;
    additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    fieldsAdd=fieldnames(additional);
    additionalMat=[];
    for ii=1:length(fieldsAdd)
        additionalMat=[additionalMat,additional.(fieldsAdd{ii})];
    end
    population.objective=objValue;
    population.additional=additionalMat;
    
    
end

function [objValue,additional]=CutCellFlow(paramsurrogate,member,loop)
    boundaryLoc=member.location;
    
    [obj]=CutCellFlow_Handler(paramsurrogate,boundaryLoc);
    [~,areaAdd]=LengthArea(member,loop);
    objValue=obj.cd;
    
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
    
end

function [objValue,additional]=LengthArea(member,loop)
    for ii=1:length(loop)
        points=loop(ii).snaxel.coord(1:end-1,:);
        [A(ii)]=abs(CalculatePolyArea(points));
        vec=points([end,1:end-1],:)-points;
        L(ii)=sum(sqrt(sum(vec.^2,2)));
        t(ii)=max(points(:,2))-min(points(:,2));
        xMin(ii)=min(points(:,1));
        xMax(ii)=max(points(:,1));
        
    end
    objValue=sum(A)/sum(L);
    
    additional.A=sum(A);
    additional.L=sum(L);
    additional.t=sum(t);
    additional.c=max(xMax)-min(xMin);
    additional.tc=additional.t/additional.c;
end

%% Execute snakes

function [snaxel,snakposition,snakSave,looprestart,restartsnake,outinfo]...
        =ExecuteSnakes_Surrogate(gridrefined,looprestart,baseGrid,connectstructinfo...
        ,param,paramspline,outinfo,nProf)
    % Executes the snakes edge detection process
    
    
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
    
    [outinfo]=SurrogateOutput('profile',[],outinfo,nProf,looprestart,restartsnake,snakSave);
    
    [textOut2,~]=evalc('PrintEnd(procStr,2,tStart)');
    fprintf([textOut1,textOut,textOut2])
end

% Print to screen functions

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

% Subdivision process

function [loop]=SubdivisionSurface_Snakes(loop,refineSteps,param,paramspline)
    % Function taking in a closed loop of vertices and applying the subdivision
    % process
    % typBOund is te type of boundary that is expected, it can either be the
    % string 'vertex' (default) or the string 'snaxel' to show that the snaxel
    % process has been used
    varExtract={'typeBound','subdivType','resampleSnak'};
    [typeBound,subdivType,resampleSnak]=ExtractVariables(varExtract,param);
    if ~exist('typeBound','var'), typeBound='vertex'; end
    
    for ii=1:length(loop)
        startPoints=loop(ii).(typeBound).coord;
        loop(ii).isccw=CCWLoop(startPoints);
        [newPoints,projPoints]=SubDivision(startPoints,refineSteps,subdivType);
        loop(ii).subdivision=newPoints;
        loop(ii).subdivspline=newPoints;
        if resampleSnak
            resampPoints=ResampleSpline(projPoints,paramspline);
            loop(ii).subdivspline=resampPoints;
        end
        %newPoints=SubSurfBSpline(startPoints(1:end-2,:),refineSteps);
        
    end
end


% Handle full results request

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


