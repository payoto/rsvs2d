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
    
    
    boundstr{1}='boundaryis0'; %'boundaryis0'
    boundstr{2}='solidnotIn0';
    boundstr{3}='0bound';
    if ~exist('caseString','var'),caseString='WeirdShape';end
    
    [passDomBounds,passGridSteps,refineSteps,passPadding...
            ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
            ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=MainInputVar(caseString);
    w = whos;
    for a = 1:length(w) 
        saveParam.(w(a).name) = eval(w(a).name); 
    end
    
    % Execution of subroutines
    
    [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(...
        passDomBounds,passGridSteps,passPadding,...
        typDat,loadLogical,isCheckRes,boundstr);
    [gridrefined,connectstructinfo,unstructuredrefined,looprefined]=...
        ExecuteGridRefinement(unstructReshape,refineGrid,boundstr,typeRefine,execTest);
    if useSnakes
        [snaxel,snakposition,snakSave,loop,cellCentredGrid]=ExecuteSnakes(unstructuredrefined,looprefined,...
            unstructured,connectstructinfo,snakesPlotInterval,snakesSteps,makeMov,boundstr);
    end
    
    
    loop=SubdivisionSurface(loop,refineSteps,typeBound);
    CheckResults(unstructured,loop,typeBound)
    BoundaryOutput(loop)
    if makeMov
        fps=5;
        quality=100;
        movStruct=[snakSave(:).movFrame];
        MakeVideo(movStruct,fps,quality,typDat);
    end
    TecplotOutput('TecPlot',typDat,unstructReshape,unstructuredrefined,snakSave,connectstructinfo)
    %OutPutBinaryResults(snakSave,saveParam,typDat)
end

%% Top Level Execution processes

function [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(...
        passDomBounds,passGridSteps,passPadding,...
        typDat,loadLogical,isCheckRes,boundstr)
    % Executes the Grid Initialisation process
    
    disp('GENERATION PROCESS START')
    t1=now;
    [unstructured,loop,unstructReshape]=...
        GridInitialisation(passDomBounds,passGridSteps,passPadding,...
        typDat,loadLogical,isCheckRes,boundstr);
    
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GENERATION PROCESS end')
    
    
end

function [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        ExecuteGridRefinement(unstructReshape,refineGrid,boundstr,typeRefine,execTest)
    % Executes the Grid Initialisation process
    
    disp('GRID REFINEMENT START')
    t1=now;
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,refineGrid,boundstr,typeRefine,execTest);
    
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GRID REFINEMENT end')
    
    
end

function [snaxel,snakposition,snakSave,loop,cellCentredGrid]=ExecuteSnakes(unstructured,loop,...
        oldGrid,connectionInfo,plotInterval,numSteps,makeMov,boundstr)
    % Executes the snakes edge detection process
    %
    
    t1=now;
    disp('SNAKE PROCESS START')
    [snaxel,snakposition,snakSave,loopsnaxel,cellCentredGrid]=Snakes(unstructured,loop,...
        oldGrid,connectionInfo,plotInterval,numSteps,makeMov,boundstr);

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

%% Output Functions to dat file

function []=BoundaryOutput(loop)
    
    % trim loops and extract data
    loopout=TrimLoops(loop);
    
    
    % format numeric data to printable strings
    cellLoops=DataToString(loopout);
    
    % print string data to file
    WriteToFile(cellLoops)
    
end

function [loopout]=TrimLoops(loop)
    % function extracting the data that must be written to the boundary.dat
    % file
    
    nLoop=length(loop);
    
    for ii=1:nLoop
        if loop(ii).isccw
            loopout.surf(ii).coord=loop(ii).subdivision(1:end-2,:);
        else
            loopout.surf(ii).coord=loop(ii).subdivision(end-2:-1:1,:);
        end
        loopout.surf(ii).nvertex=length(loopout.surf(ii).coord(:,1));
        loopout.surf(ii).nfaces=loopout.surf(ii).nvertex;
    end
    loopout.total.nvertex=sum([loopout.surf(:).nvertex]);
    loopout.total.nfaces=sum([loopout.surf(:).nfaces]);
    loopout.total.nloops=nLoop;
    
end

function cellLoops=DataToString(loopout)
    % Transforms the data organised in loopout into a cell array of strings
    % ready to be written to a file
    
    cellLoops{1}=int2str(loopout.total.nloops);
    cellLoops{2}=int2str([loopout.total.nvertex, loopout.total.nfaces]);
    kk=3;
    
    % run through the different surfaces
    for ii=1:loopout.total.nloops
        cellLoops{kk}=int2str(loopout.surf(ii).nvertex);
        kk=kk+1;
        % Write an array of numbers in the correct format
        for jj=1:length(loopout.surf(ii).coord(:,1))
            cellLoops{kk}='';
            for ll=1:length(loopout.surf(ii).coord(jj,:))
                cellLoops{kk}=[cellLoops{kk},...
                    num2str(loopout.surf(ii).coord(jj,ll),8),'  '];
            end
            kk=kk+1;
        end
        
    end
    
end

function []=WriteToFile(cellLoops)
    % writes the data to the file
    [FID]=OpenValidFile();
    
    for ii=1:length(cellLoops)
        for jj=1:length(cellLoops{ii}(:,1))
            fprintf(FID,[cellLoops{ii}(jj,:),'\n']);
        end
    end
    
    
    fclose(FID);
end

function [FID]=OpenValidFile(writeDirectory)
    % Creates a file in the current directory to write data to.
    
    datestr(now)
    if ~exist('writeDirectory','var')
        cdDir=cd;
        backSlash=regexp(cdDir,'\');
        writeDirectory=[cdDir,'boundaries'];
    end
    system(['md "',writeDirectory,'"']);
    
    fileName=['boundary_',datestr(now,30),'.dat'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

%% General Utility Functions

function []=MakeVideo(movStruct,fps,quality,datType)
    
    fileName=['Vid_',datType,'_',datestr(now,30),'.avi'];
    stampedSubFolders=['VideoArchive_',datestr(now,'yyyy_mm')];
    targetDir=[cd,'\Videos\',stampedSubFolders,'\'];
    system(['md "',targetDir,'"']);
    writerObj = VideoWriter([targetDir,fileName]);
    writerObj.FrameRate=fps;
    writerObj.Quality=quality;
    open(writerObj)
    writeVideo(writerObj,movStruct)
    close(writerObj)
end

function []=OutPutBinaryResults(snakSave,paramStruct,datType,optionalSubFolder)
    if ~exist('optionalSubFolder','var'),optionalSubFolder='';end
    
    fileName=['Variables_',datType,'_',datestr(now,30),'.mat'];
    stampedSubFolders=['DataArchive_',datestr(now,'yyyy_mm')];
    targetDir=[cd,'\Results\',optionalSubFolder,'\',stampedSubFolders,'\'];
    system(['md "',targetDir,'"']);
    
    save([targetDir,fileName],'snakSave','paramStruct');
    
end

function []=template()
    
end
