%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [unstructured,loop,unstructReshape,snakSave,param]=Main(caseString,restart)
    % Main function for the execution of the Subdivision process
    
    close all
    clc
    
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    include_PostProcessing
    include_Mex_Wrapper
    
    diaryFile=[cd,'\Result_Template\Latest_Diary.log'];
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    
    
    if ~exist('caseString','var'),caseString='WeirdShape';end
    if ~exist('restart','var'),restart=false;end
    
    if ~restart
        [param,unstructured,unstructuredrefined,loop,connectstructinfo...
            ,snakSave,unstructReshape]=StandardRun(caseString);
    else
        [param,unstructured,unstructuredrefined,loop,connectstructinfo...
            ,snakSave,unstructReshape]=RestartRun(caseString);
    end
    varExtract={'useSnakes','typeBound','refineSteps'};
    [useSnakes,typeBound,refineSteps]=ExtractVariables(varExtract,param);
    
    % Restart structure generation
    restartstruct.unstructuredrefined=unstructuredrefined;
    restartstruct.param=param;
    restartstruct.connectstructinfo=connectstructinfo;
    restartstruct.unstructured=unstructured;
    if useSnakes
        restartstruct.snakrestart.snaxel=snakSave(end).snaxel;
        restartstruct.snakrestart.cellCentredGid=snakSave(end).cellCentredGrid;
        restartstruct.snakrestart.insideContourInfo=snakSave(end).insideContourInfo;
    end
    % Post processes
    loop=SubdivisionSurface(loop,refineSteps,typeBound);
    CheckResults(unstructured,loop,typeBound)
    
    tecoutstruct.baseGrid=unstructReshape;
    tecoutstruct.fineGrid=unstructuredrefined;
    tecoutstruct.snakSave=snakSave;
    tecoutstruct.connectstructinfo=connectstructinfo;
    
    diary off
    ManageOutputResults(param,loop,tecoutstruct,restartstruct);
    %TecplotOutput(unstructReshape,unstructuredrefined,snakSave,connectstructinfo)
    %OutPutBinaryResults(snakSave,saveParam,typDat)
end

%% Top Level Execution processes

function [param,unstructured,unstructuredrefined,loop,connectstructinfo...
        ,snakSave,unstructReshape]=StandardRun(caseString)
    
    [param]=structInputVar(caseString);
    param.general.restart=false;
    varExtract={'useSnakes'};
    [useSnakes]=ExtractVariables(varExtract,param);
    
    
    % Execution of subroutines
    
    [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(param);
    
    [gridrefined,connectstructinfo,unstructuredrefined,looprefined]=...
        ExecuteGridRefinement(unstructReshape,param);
    snakSave=[];
    if useSnakes
        [snaxel,snakposition,snakSave,loop,cellCentredGrid]=ExecuteSnakes(unstructuredrefined,looprefined,...
            unstructured,connectstructinfo,param);
    end
    
end

function [param,unstructured,unstructuredrefined,loop,connectstructinfo...
        ,snakSave,unstructReshape]=RestartRun(caseStr)
    
    load([caseStr,'.mat'])
    param.general.restart=true;
    
    FID=fopen([cd,'\param_Restart.m'],'w');
    GenerateParameterFile(FID,param,now,'Restart')
    
    disp(' ')
    disp('Edit Parameter File in Notepad++, press Any key in the MATLAB command window')
    disp(' once it is saved to resume execution.')
    system(['"C:\Program Files (x86)\Notepad++\notepad++.exe"  ',...
        '"',[cd,'\param_Restart.m'],'"']);
    pause
    
    param_Restart
    varExtract={'useSnakes'};
    [useSnakes]=ExtractVariables(varExtract,param);
    if useSnakes
        [snaxel,snakposition,snakSave,loop,cellCentredGrid]=ExecuteSnakes(unstructuredrefined,snakrestart,...
            unstructured,connectstructinfo,param);
    end
    unstructReshape=ModifUnstructured(unstructured);
end

function [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(param)
    % Executes the Grid Initialisation process
    
    disp('GENERATION PROCESS START')
    t1=now;
    [unstructured,loop,unstructReshape]=...
        GridInitialisationV2(param);
    
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GENERATION PROCESS end')
    
    
end

function [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        ExecuteGridRefinement(unstructReshape,param)
    % Executes the Grid Initialisation process
    
    disp('GRID REFINEMENT START')
    t1=now;
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GRID REFINEMENT end')
    
    
end

function [snaxel,snakposition,snakSave,loop,cellCentredGrid]=ExecuteSnakes(unstructured,loop,...
        oldGrid,connectionInfo,param)
    % Executes the snakes edge detection process
    %
    
    t1=now;
    disp('SNAKE PROCESS START')
    [snaxel,snakposition,snakSave,loopsnaxel,cellCentredGrid]=Snakes(unstructured,loop,...
        oldGrid,connectionInfo,param);

    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('SNAKE PROCESS END')

    if length(loopsnaxel)==length(loop)
        for ii=1:length(loopsnaxel)
            loop(ii).snaxel=loopsnaxel(ii).snaxel;
        end
    else
        loop=loopsnaxel;
    end
end


%% Subdivision process

function [loop]=SubdivisionSurface(loop,refineSteps,typeBound)
    % Function taking in a closed loop of vertices and applying the subdivision
    % process
    % typBOund is te type of boundary that is expected, it can either be the
    % string 'vertex' (default) or the string 'snaxel' to show that the snaxel
    % process has been used
    
    if ~exist('typeBound','var'), typeBound='vertex'; end
    
    for ii=1:length(loop)
        startPoints=loop(ii).(typeBound).coord;
        loop(ii).isccw=CCWLoop(startPoints);
        newPoints=SubSurfChainkin(startPoints(1:end-2,:),refineSteps);
        loop(ii).subdivision=newPoints;
    end
end

function [newPoints]=SubSurfChainkin2(startPoints,refineSteps)
    % Implements a Chainkin subdivision process

    
    chainkinMask=[0.25,0.75,0.75,0.25]';
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        subMask=zeros((numPoints*2+2),numPoints);
        
        % create the right sized mask
        for ii=1:numPoints
            jj=(1+(ii-1)*2);
            subMask(jj:jj+3,ii)=chainkinMask;
        end
        
        newPoints=subMask*startPoints;
        newPoints(1:2,:)=[];
        newPoints(end-1:end,:)=[];
        startPoints=newPoints;
    end
end

function [newPoints]=SubSurfChainkin(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    TEisLeft=0;
    chainkinNoCorn=zeros([4,3]);
    chainkinNoCorn(1:4,2)=[0.25;0.75;0.75;0.25];
    
    chainkinCorn=zeros([5,3]);
    chainkinCorn(1:5,1)=[0;0.25;0;0;0];
    chainkinCorn(1:5,2)=[0.25;0.5;1;0.5;0.25];
    chainkinCorn(1:5,3)=[0;0;0;0.25;0];
    
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        isCorner=DetectTrailingEdge(startPoints,TEisLeft);
        cumCorner=cumsum(isCorner);
        numNewPoints=(numPoints*2+cumCorner(end));
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii-1;
            jStart=ii*2+cumCorner(ii+1)-isCorner(ii+1);
            
            if isCorner(ii+1)
                nJ=5;
                chainkinMask=chainkinCorn;
            else
                nJ=4;
                chainkinMask=chainkinNoCorn;
            end
            indX=zeros(1,3);
            indY=zeros(1,nJ);
            for iLoop=1:3
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=chainkinMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPoints=newPoints;
        
    end
end
function isCorner=DetectTrailingEdge(coord,TEisLeft)
    TEisLeft=(TEisLeft-0.5)*2;
    testLocMin=((TEisLeft*(coord([2:end,1],1)-coord(:,1)))>0) ...
        & ((TEisLeft*(coord([end,1:end-1],1)-coord(:,1)))>0);
    isCorner=testLocMin;
end
%% Plot Functions
function []=CheckResults(unstructured,loop,typeBound)
    global nDim domainBounds
    
    if nDim==2
        figh=figure;
        axh=axes;
        hold on
        
        colString='bgcmyk';
        
        isEdgeIndex=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeIndex)
            PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'bo')
        end
        
        isEdgeIndex=find(unstructured.edge.boundaryis0);
        for ii=1:length(isEdgeIndex)
            PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'b-')
        end
        
        
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        
        for ii=1:length(loop)
            [~,colIndex]=IntegerQuotient(ii,length(colString));
            colIndex=colIndex+1;
            PlotLoop(figh,axh,loop,ii,[colString(colIndex),'--'],typeBound)
            PlotSubDiv(figh,axh,loop,ii,[colString(colIndex),'-'])
        end
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    end
    
end

function []=PlotEdge(figh,axh,unstructured,indexEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(indexEdge,:);
    coord=unstructured.vertex.coord(vertices,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotLoop(figh,axh,loop,indexLoop,format,typeBound)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).(typeBound).coord;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotSubDiv(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).subdivision;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotCell(figh,axh,unstructured,indexCell,format)
    figure(figh)
    axes(axh)
    
    
    coord=unstructured.cell.coord(indexCell,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end


%% General Utility Functions

function []=OutPutBinaryResults(snakSave,paramStruct,datType,optionalSubFolder)
    if ~exist('optionalSubFolder','var'),optionalSubFolder='';end
    
    fileName=['Variables_',datType,'_',datestr(now,30),'.mat'];
    stampedSubFolders=['DataArchive_',datestr(now,'yyyy_mm')];
    targetDir=[cd,'\..\results\',optionalSubFolder,'\',stampedSubFolders,'\'];
    system(['md "',targetDir,'"']);
    
    save([targetDir,fileName],'snakSave','paramStruct');
    
end

function []=template()
    
end
