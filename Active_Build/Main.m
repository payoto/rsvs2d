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


function [unstructured,loop,unstructReshape,snakSave]=Main(caseString)
    % Main function for the execution of the Subdivision process
    
    close all
    clc
    
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    
    boundstr{1}='boundaryis0'; %'boundaryis0'
    boundstr{2}='solidnotIn0';
    boundstr{3}='0bound';
    if ~exist('caseString','var'),caseString='WeirdShape';end
    
    [param]=structInputVar(caseString);
    
    varExtract={'useSnakes','typeBound','refineSteps','makeMov','typDat'};
    [useSnakes,typeBound,refineSteps,makeMov,typDat]=ExtractVariables(varExtract,param);
    
    
    w = whos;
    for a = 1:length(w) 
        saveParam.(w(a).name) = eval(w(a).name); 
    end
    
    % Execution of subroutines
    
    [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(param);
    
    [gridrefined,connectstructinfo,unstructuredrefined,looprefined]=...
        ExecuteGridRefinement(unstructReshape,param);
    
    if useSnakes
        [snaxel,snakposition,snakSave,loop,cellCentredGrid]=ExecuteSnakes(unstructuredrefined,looprefined,...
            unstructured,connectstructinfo,param);
    end
    
    
    loop=SubdivisionSurface(loop,refineSteps,typeBound);
    CheckResults(unstructured,loop,typeBound)
    
    tecoutstruct.baseGrid=unstructReshape;
    tecoutstruct.fineGrid=unstructuredrefined;
    tecoutstruct.snakSave=snakSave;
    tecoutstruct.connectstructinfo=connectstructinfo;
    
    ManageOutputResults(param,loop,tecoutstruct)
    %TecplotOutput(unstructReshape,unstructuredrefined,snakSave,connectstructinfo)
    %OutPutBinaryResults(snakSave,saveParam,typDat)
end

%% Top Level Execution processes

function [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(param)
    % Executes the Grid Initialisation process
    
    disp('GENERATION PROCESS START')
    t1=now;
    [unstructured,loop,unstructReshape]=...
        GridInitialisation(param);
    
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
        newPoints=SubSurfChainkin(startPoints,refineSteps);
        loop(ii).subdivision=newPoints;
    end
end

function [newPoints]=SubSurfChainkin(startPoints,refineSteps)
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
