%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%                - Snakes -
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main execution functions
function [snaxel,snakposition,snakSave,loopsnaxel,cellCentredGrid]=Snakes(refinedGriduns,loop,...
        oldGridUns,connectionInfo,plotInterval,numSteps,makeMov,boundstr)
    
    global arrivalTolerance
    arrivalTolerance=0.01;
    
    global unstructglobal
    unstructglobal=refinedGriduns;
    
    global maxStep
    maxStep=0.5;
    maxDt=0.9;
    
    forceparam.maxForceVel=1;
    forceparam.bendingVelInfluence=0.5;
    forceparam.tensVelInfluence=1;
    forceparam.maxVelRatio=1;
    
    if ~exist('numSteps','var'),numSteps=50;end
    if ~exist('plotInterval','var'),plotInterval=ceil(numSteps/10);end
    [refinedGrid]=ModifUnstructured(refinedGriduns);
    [oldGrid]=ModifUnstructured(oldGridUns);
    
    debugPlot=[53:58];%[1:10:400,2, 4,6,8];
    mergeTopo=true;
    
    [snaxel,snakposition,snakSave,loopsnaxel,cellCentredGrid]=...
        RunSnakesProcess(refinedGriduns,refinedGrid,loop,...
        oldGrid,oldGridUns,connectionInfo,plotInterval,numSteps,debugPlot,...
        arrivalTolerance,maxStep,maxDt,mergeTopo,makeMov,boundstr,forceparam);
    figure,semilogy(1:length(snakSave),[snakSave(:).currentConvVolume])
    title('Volume error')
    ylabel('Root Mean squared error on volume convergence')
    xlabel('number of iterations')
end

function [snaxel,snakposition,snakSave,loopsnaxel,cellCentredGrid]=...
        RunSnakesProcess(refinedGriduns,refinedGrid,loop,...
        oldGrid,oldGridUns,connectionInfo,plotInterval,numSteps,debugPlot,...
        arrivalTolerance1,maxStep,maxDt,mergeTopo,makeMov,boundstr,forceparam)
    % Main execution container for Snakes
    
    global arrivalTolerance
    arrivalTolerance=arrivalTolerance1;
    dtMin=1;
    decayCoeff=0.15;
    convLevel=10^-8;
    
    
    disp(['    Start initialisation'])
        tStepStart=now;
    [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        StartSnakeProcess(refinedGriduns,refinedGrid,loop,...
        oldGrid,connectionInfo,mergeTopo,boundstr);
    tStepEnd=now;
    disp(['   Initialisation time:',datestr(tStepEnd-tStepStart,'HH:MM:SS:FFF')]);
        
    tStart=now;
    movFrame=struct([]);
    for ii=1:numSteps
        disp(' ')
        disp(['Start step ',num2str(ii)])
        tStepStart=now;
        
        %arrivalTolerance=arrivalTolerance1*exp(-decayCoeff*ii);
        % snaxel properties calculation
        [snakposition]=PositionSnakes(snaxel,refinedGriduns);
        [snakposition]=SnaxelNormal2(snaxel,snakposition);
        [volumefraction,coeffstructure,cellCentredGridSnax]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
            cellCentredGrid,insideContourInfo);
        [snaxel,snakposition,snaxelmodvel]=VelocityCalculationVolumeFraction(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        % visualise results
        
        if ((round(ii/plotInterval)==ii/plotInterval) && plotInterval) || sum(ii==debugPlot)
            [movFrame]=CheckResults(ii,refinedGriduns,oldGrid,snakposition,snaxelmodvel,makeMov,volumefraction);
            [movFrame]=CheckResults(ii,refinedGriduns,oldGrid,snakposition,snaxel,makeMov,volumefraction);
        
        end
        
        [convergenceCondition,currentConvVelocity,currentConvVolume]=...
            ConvergenceTest(snaxel,volumefraction,convLevel);
        if convergenceCondition
            break
        end
        
        [dt,dtSnax,maxDist]=TimeStepCalculation(snaxel,maxStep,maxDt,dtMin);
        
        % Save and exit conditions
        snakSave(ii).snaxel=snaxel;
        snakSave(ii).dt=dt;
        snakSave(ii).snakposition=snakposition;
        snakSave(ii).volumefraction=volumefraction;
        snakSave(ii).cellCentredGrid=cellCentredGridSnax;
        snakSave(ii).currentConvVelocity=currentConvVelocity;
        snakSave(ii).currentConvVolume=currentConvVolume;
        snakSave(ii).movFrame=movFrame;
        
        snaxel=SnaxelDistanceUpdate(snaxel,dt,dtSnax,maxDist);
        
        [snakposition]=PositionSnakes(snaxel,refinedGriduns);
        [snakposition]=SnaxelNormal2(snaxel,snakposition);
        [volumefraction,coeffstructure,cellCentredGridSnax]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
            cellCentredGrid,insideContourInfo);
        [snaxel,snakposition,snaxelmodvel]=VelocityCalculationVolumeFraction(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        
        [snaxel]=FreezingFunction(snaxel,borderVertices,mergeTopo);
        % Snaxel Repopulation In both directions
        suba=ArrivalCondition(snaxel);
        [savTest(ii).finishedSnakes]=[snaxel(suba).index];
        [snaxel,insideContourInfo]=SnaxelRepopulate(refinedGriduns,snaxel,...
            insideContourInfo);
        [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
        [snaxelrev,insideContourInfoRev]=ReverseSnaxelInformation(snaxel,...
            insideContourInfo,refinedGriduns);
        
        
        suba=ArrivalCondition(snaxelrev);
        [savTest(ii).finishedSnakesrev]=[snaxelrev(suba).index];
        [snaxelrev,insideContourInfoRev]=SnaxelRepopulate(refinedGriduns,snaxelrev,...
            insideContourInfoRev);
        [snaxelrev,insideContourInfoRev]=SnaxelCleaningProcess(snaxelrev,insideContourInfoRev);
        [snaxel,insideContourInfo]=ReverseSnaxelInformation(snaxelrev,...
            insideContourInfoRev,refinedGriduns);
        
        
        % Topology Trimming, Merging and Freezing
        [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
        [snakposition]=PositionSnakes(snaxel,refinedGriduns);
        [snaxel]=FreezingFunction(snaxel,borderVertices,mergeTopo);
        [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo);
        
        
        
        
        nSnax=length(snaxel);
        frozen=sum([snaxel(:).isfreeze]);
        tStepEnd=now;
        disp(['   Step time:',datestr(tStepEnd-tStepStart,'HH:MM:SS:FFF')]);
        
        if nSnax==frozen
            break
        end
        
    end
    tEnd=now;
    disp(['Iteration Time:',datestr(tEnd-tStart,'HH:MM:SS:FFF')]);
    disp(['Volume converged to ',num2str(currentConvVolume,'%.5e')])
    [snaxel,snakposition,loopsnaxel]=FinishSnakes(snaxel,...
        borderVertices,refinedGriduns);
    
    %[snakSave]=ReformatSnakSave(snakSave);
end

function [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        StartSnakeProcess(refinedGriduns,refinedGrid,loop,...
        oldGrid,connectionInfo,mergeTopo,boundstr)
    
    [cellCentredGrid]=CellCentredGridInformation(refinedGrid);
    [volfracconnec]=VolumeFractionConnectivity(oldGrid,...
        connectionInfo,cellCentredGrid,refinedGrid);
    
    
    insideContourInfo=refinedGriduns.edge.(boundstr{2});
    disp('Find Border Vertices')
    [borderVertices]=FindBorderVertex(refinedGriduns);
    disp('Initialise Snaxel Grid')
    
    % Inside contour info will change depending on the type of contour under
    % consideration (0 contour or 1 contour)
    [snaxel,insideContourInfo]=SnaxelInitialisation(refinedGriduns,loop,insideContourInfo,boundstr);
    
    [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo);
    [snaxel]=FreezingFunction(snaxel,borderVertices,mergeTopo);
end

function [snaxel,snakposition,loopsnaxel]=FinishSnakes(snaxel,...
        borderVertices,refinedGriduns)
    
    disp('Finished Iterations , starting Post Process')
    [snaxel]=FreezingFunction(snaxel,borderVertices);
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    CheckResultsLight(refinedGriduns,snakposition,snaxel)
    disp('Creating Snaxel Loops')
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel);
    
end

function [snakSave]=ReformatSnakSave(snakSave)
    % Reformats the normal vectors in snak Save
    
    for ii=1:length(snakSave)
        for jj=1:length(snakSave(ii).snakposition)
            snakSave(ii).snakposition(jj).normvector=...
                [snakSave(ii).snakposition(jj).normvector{1};...
                snakSave(ii).snakposition(jj).normvector{2}];
        end
        
    end
end

function [snaxelRev,insideContourInfoRev]=ReverseSnaxelInformation(snaxel,...
        insideContourInfo,unstructured)
    % Reverses all the snaxel information
    [snaxelRev]=ReverseSnakes(snaxel);
    [insideContourInfoRev]=ReverseInsideContourInfo(snaxel,...
        insideContourInfo,unstructured);
    
end

function [convergenceCondition,currentConvVelocity,currentConvVolume]=...
        ConvergenceTest(snaxel,volumefraction,convLevel)
    % Calculates wether the snaxel population has converged
    
    vSnax=[snaxel(:).v];
    currentConvVelocity=(sqrt(sum(vSnax.^2))/length(vSnax));
    conditionVelocity=currentConvVelocity<convLevel;
    
    diffVolFrac=[volumefraction(:).targetfill]-[volumefraction(:).volumefraction];
    currentConvVolume=(sqrt(sum(diffVolFrac.^2))/length(diffVolFrac));
    conditionVolume=currentConvVolume<convLevel;
    
    %convergenceCondition= conditionVolume && conditionVelocity;
    convergenceCondition= conditionVolume;
end

%% Snaxel Initialisation

function [snaxel,insideContourInfo]=SnaxelInitialisation(unstructured,loop,insideContourInfo,boundstr)
    
    for ii=1:length(loop)
        if ii==1
            snaxelIndexStart=0;
        else
            snaxelIndexStart=max([snaxel(:).index]);
        end
        [loopsnakes(ii).snaxels,insideContourInfo]=SnaxelLoop(unstructured,loop(ii),...
            snaxelIndexStart,insideContourInfo,boundstr);
        if ii==1
            snaxel=loopsnakes(ii).snaxels;
        else
            snaxel=[snaxel,loopsnakes(ii).snaxels];
        end
    end
    
end

function [snaxel,insideContourInfo]=SnaxelLoop(unstructured,loop,...
        snaxelIndexStart,insideContourInfo,boundstr)
    
    loopVertInd=loop.vertex.index(1:end-2); %RemoveIdenticalEntries(loop.vertex.index);
    loopVertSub=FindObjNum([],loopVertInd,unstructured.vertex.index);
    vertCoord=loop.vertex.coord(1:end-2,:);
    isLoopCCW=CCWLoop(vertCoord);
    if isLoopCCW
        initVertexIndex=loopVertInd;
    else
        %initVertexIndex=loop.vertex.index(end-2:-1:1);
        initVertexIndex=loopVertInd(end:-1:1);
    end
    %initVertexIndex=loop.vertex.index(1:end-2);
    %initVertexIndex=loop.vertex.index([2:3,5:end-3]);
    edgeVertIndex=unstructured.edge.vertexindex;
    edgeIndex=unstructured.edge.index;
    loopEdgeIndex=loop.edge.index;
    
    [snaxel,~]=InitialSnaxelStructure(initVertexIndex,edgeVertIndex,...
        edgeIndex,loopEdgeIndex,snaxelIndexStart);
    
    [delIndex]=FindInsideSnaxels(snaxel,insideContourInfo);
    switch boundstr{3}
        case '1bound'
            loopEdgeSubs=FindObjNum([],loopEdgeIndex,edgeIndex);
        case '0bound'
            snaxInd=[snaxel(:).index];
            keepIndex=delIndex;
            keepSub=FindObjNum([],keepIndex,snaxInd);
            logDelInd=true(size(snaxInd));
            logDelInd(keepSub)=false;
            delIndex=snaxInd(logDelInd);
            loopEdgeSubs=[];
            snaxel=ReverseSnakes(snaxel);
        otherwise
            error('Invalid boundstr flag')
            
    end
    
    snaxel=DeleteSnaxel(snaxel,delIndex);
    
    [snaxel]=TestSnaxelLoopDirection(snaxel);
    
    insideContourInfo(loopEdgeSubs)=1;
    
end

function [snaxel]=TestSnaxelLoopDirection(snaxel)
    
    global unstructglobal
    snakpos=PositionSnakes(snaxel,unstructglobal);
    coord=vertcat(snakpos(:).coord);
    [leftMostCorner]=LeftMostCorner(coord);
    leftCornerDirection=snakpos(leftMostCorner).vector;
    testAngles=ExtractAngle360(leftCornerDirection,[1 0]);
    
    nextSnaxInd=snaxel(leftMostCorner).snaxnext;
    precSnaxInd=snaxel(leftMostCorner).snaxprec;
    nextSnaxSub=FindObjNum(snaxel,nextSnaxInd);
    precSnaxSub=FindObjNum(snaxel,precSnaxInd);

    precVec=coord(precSnaxSub,:)-coord(leftMostCorner,:);
    nextVec=coord(nextSnaxSub,:)-coord(leftMostCorner,:);
    precAngle=ExtractAngle360(precVec,[-1 -1]);
    nextAngle=ExtractAngle360(nextVec,[-1 -1]);


    if precAngle<nextAngle
        isCCW=true;
    elseif precAngle>nextAngle
        isCCW=false;
    else
        isCCW=[];
    end
    
    if sum(testAngles>pi/2)==0
        if isCCW
            [snaxel]=ReverseSnakesConnection(snaxel);
        end
    else
        if ~isCCW
            [snaxel]=ReverseSnakesConnection(snaxel);
        end
    end
    
end

function [snaxel,cellSimVertex]=InitialSnaxelStructure(initVertexIndex,edgeVertIndex,...
        edgeIndex,invalidEdgeIndex,snaxelIndexStart,connectivity)
    % Creates the basic snaxel structure, no removal of unworthy snaxels is
    % performed
    
    if ~exist('connectivity','var'); isConnec=false; else isConnec=true; end
    [mIVI,~]=size(initVertexIndex);
    
    kk=0;
    kk2=0;
    cellSimVertex{mIVI}=[];
    snaxel=[];
    
    
    for ii=1:mIVI
        vertexEdges=FindEdgesIndex(initVertexIndex(ii),edgeVertIndex,edgeIndex);
        snaxelEdges=FindSnaxelEdges(vertexEdges,invalidEdgeIndex);
        if numel(snaxelEdges)>0
            [kk,cellSimVertex{ii},snaxelNotOrdered]=GenerateVertexSnaxel(snaxelEdges,kk,...
                snaxelIndexStart,initVertexIndex(ii), edgeVertIndex,edgeIndex);
            
            % Now Fine
            if numel(invalidEdgeIndex)>1
                baseEdgeExploit=zeros(length(invalidEdgeIndex),1);
                for jj=1:length(invalidEdgeIndex)
                    if sum(vertexEdges==invalidEdgeIndex(jj))
                        baseEdgeExploit(jj)=1;
                    end
                end
                baseEdgeExploitTest=baseEdgeExploit+[baseEdgeExploit(2:end);baseEdgeExploit(1)];
                baseEdgeExploit=find(baseEdgeExploitTest==2);
            else
                baseEdgeExploit=1;
            end
            
            connecOrder=CCWOrderAroundNode(snaxelNotOrdered,invalidEdgeIndex(baseEdgeExploit));
            generateOrder=FindObjNum(snaxelNotOrdered,connecOrder);
            
            [kk2,cellSimVertex{ii},newsnaxel]=GenerateVertexSnaxel...
                (snaxelEdges(generateOrder),kk2,snaxelIndexStart,initVertexIndex(ii),...
                edgeVertIndex,edgeIndex);
            snaxel=[snaxel,newsnaxel];
        end
    end
    if isConnec
        snaxel(1).connectivity(1)=connectivity(1);
        snaxel(1).snaxprec=connectivity(1);
        
        snaxel(end).connectivity(2)=connectivity(2);
        snaxel(end).snaxnext=connectivity(2);
    else
        snaxel(1).connectivity(1)=snaxel(end).index;
        snaxel(1).snaxprec=snaxel(end).index;
        snaxel(end).connectivity(2)=snaxel(1).index;
        snaxel(end).snaxnext=snaxel(1).index;
    end
end

function [kk,cellSimVertex,snaxel]=GenerateVertexSnaxel(snaxelEdges,kk,...
        snaxelIndexStart,initVertexIndexSingle, edgeVertIndex,edgeIndex)
    
    numSE=length(snaxelEdges); % provides information about snaxel from same vertex
    snaxelEdgesSub=FindObjNum([],snaxelEdges,edgeIndex);
    kkLocal=0;
    for jj=1:numSE
        kk=kk+1;
        kkLocal=kkLocal+1;
        snaxIndex=kk+snaxelIndexStart;
        cellSimVertex(jj)=snaxIndex;
        dist=0; % Snaxel initialisation, it starts at the vertex
        velocity=1; % Initialisation velocity
        vertexOrig=initVertexIndexSingle;
        vertexDest=edgeVertIndex(snaxelEdgesSub(jj),:);
        vertexDest(vertexDest==vertexOrig)=[];
        currEdge=snaxelEdges(jj);
        if isempty(vertexDest)
            vertexDest=0;
        end
        snaxPrec=snaxIndex-1;
        snaxNext=snaxIndex+1;
        snaxel(kkLocal)=SnaxelStructure(snaxIndex,dist,velocity,vertexOrig,...
            vertexDest,snaxPrec,snaxNext,currEdge);
    end
    
end

function edges=FindEdgesIndex(vertexIndex,edgeVertIndex,edgeIndex)
    % finds the edges connected to a vertex
    
    edgesPosition=find(sum(edgeVertIndex==vertexIndex,2));
    edges=edgeIndex(edgesPosition);
end

function snaxelEdges=FindSnaxelEdges(vertexEdges,invalidEdgeIndex)
    % finds the edges connected to a vertex which will create a snaxel
    
    
    numVE=length(vertexEdges);
    isSnaxelEdge=false([numVE,1]);
    for jj=1:numVE
        isSnaxelEdge(jj)=sum(invalidEdgeIndex==vertexEdges(jj))==0;
    end
    snaxelEdges=vertexEdges(isSnaxelEdge);
end

function [snaxel]=SnaxelStructure(index,dist,velocity,vertexOrig,...
        vertexDest,snaxPrec,snaxNext,edge)
    
    snaxel.index=index;
    snaxel.d=dist;
    snaxel.v=velocity;
    snaxel.fromvertex=vertexOrig;
    snaxel.tovertex=vertexDest;
    snaxel.edge=edge;
    snaxel.connectivity=[snaxPrec snaxNext];
    snaxel.snaxprec=snaxPrec;
    snaxel.snaxnext=snaxNext;
    snaxel.isfreeze=0;
end


%% Visualisation functions

function [movFrame]=CheckResults(iter,unstructured,oldGrid,snakposition,snaxel,makeMovie,volumefraction,borderVertices)
    global nDim domainBounds
    movFrame=[];
    if nDim==2
        figh=figure('Position',[100 100 1000 900]);
        axh=axes;
        hold on
        title(['Iteration  ',int2str(iter)],'fontsize',16)
        colString='bgcmyk';
        
        isEdgeSub=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeSub)
            PlotEdge(figh,axh,unstructured,isEdgeSub(ii),'bo')
        end
        
        isEdgeSub=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeSub)
            PlotEdge(figh,axh,unstructured,isEdgeSub(ii),'b-')
        end
        if exist('borderVertices','var')
            subVert=FindObjNum([],borderVertices,unstructured.vertex.index);
            for ii=1:length(subVert)
                PlotVert(figh,axh,unstructured,subVert(ii),'ro')
            end
        end
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        PlotSnaxel(figh,axh,snakposition,snaxel)
        %PlotSnaxelLoop(figh,axh,snakposition,snaxel)
        PlotSnaxelLoopDir(figh,axh,snakposition,snaxel)
        PlotSnaxelIndex(figh,axh,snakposition)
        
        if exist('volumefraction','var')
            oldCellIndUnstructInd=[oldGrid.cell(:).index];
            
            oldCellIndUnstructSub=FindObjNum(oldGrid.cell,oldCellIndUnstructInd);
            oldCellIndVolFracSub=FindObjNum(volumefraction,...
                oldCellIndUnstructInd,[volumefraction(:).oldCellInd]);
            for ii=1:length(oldCellIndUnstructInd)
                
                coord=oldGrid.cell(oldCellIndUnstructSub(ii)).coord;
                frac=volumefraction(oldCellIndVolFracSub(ii)).volumefraction...
                    -volumefraction(oldCellIndVolFracSub(ii)).targetfill;
                %frac=volumefraction(oldCellIndVolFracSub(ii)).oldCellInd;
                PlotVolFrac(figh,axh,coord,frac)
            end
        end
        
%         [normalcontourvec]=ContourNormal2(snaxel,snakposition);
%         PlotContVec(figh,axh,snakposition,normalcontourvec)
        
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
        if makeMovie
            movFrame=getframe(figh);
        end
    end
    
end

function []=CheckResultsLight(unstructured,snakposition,snaxel)
    global nDim domainBounds
    
    if nDim==2
        figh=figure;
        axh=axes;
        hold on
        
        colString='bgcmyk';
        
        isEdgeIndex=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeIndex)
            %PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'bo')
        end
        
        isEdgeIndex=find(unstructured.edge.boundaryis0);
        for ii=1:length(isEdgeIndex)
            %PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'b-')
        end
        
        
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        %PlotSnaxel(figh,axh,snakposition)
        PlotSnaxelLoop(figh,axh,snakposition,snaxel)
        %PlotSnaxelIndex(figh,axh,snakposition)
        
        %[normalcontourvec]=ContourNormal(snaxel,snakposition);
        %PlotContVec(figh,axh,snakposition,normalcontourvec)
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    end
    
end

function []=PlotEdge(figh,axh,unstructured,subEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(subEdge,:);
    vertsub(1)=find(unstructured.vertex.index==vertices(1));
    vertsub(2)=find(unstructured.vertex.index==vertices(2));
    coord=unstructured.vertex.coord(vertsub,:);
    
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotVert(figh,axh,unstructured,subVert,format)
    figure(figh)
    %axes(axh)
    
    coord=unstructured.vertex.coord(subVert,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotSnaxel(figh,axh,snakposition,snaxel)
    % Plots the snaxels as arrows on the plot
    for ii=1:length(snakposition)
        X(ii)=snakposition(ii).coord(1);
        Y(ii)=snakposition(ii).coord(2);
        U(ii)=snakposition(ii).vector(1)*snaxel(ii).v/40;
        V(ii)=snakposition(ii).vector(2)*snaxel(ii).v/40;
    end
    figure(figh)
    axes(axh)
    quiver(X,Y,U,V,0,'r-')
    
end

function []=PlotContVec(figh,axh,snakposition,normalContVec)
    % Plots the snaxels as arrows on the plot
    snaxIndex=[snakposition(:).index];
    for ii=1:length(normalContVec)
        for jj=1:2
            workInd=normalContVec(ii).(['index',int2str(jj)]);
            %workInd=normalContVec(ii).vertex(jj);
            workSub(jj)=FindObjNum(snakposition,workInd,snaxIndex);
            coord(jj,1:2)=snakposition(workSub(jj)).coord;
        end
        coord=mean(coord);
        X(ii)=coord(1);
        Y(ii)=coord(2);
        U(ii)=normalContVec(ii).vector(1)/20;
        V(ii)=normalContVec(ii).vector(2)/20;
    end
    figure(figh)
    axes(axh)
    quiver(X,Y,U,V,0,'r-')
    
end

function []=PlotSnaxelIndex(figh,axh,snakposition)
    % Plots the snaxels as arrows on the plot
    figure(figh)
    axes(axh)
    for ii=1:length(snakposition)
        X(ii)=snakposition(ii).coord(1);
        Y(ii)=snakposition(ii).coord(2);
        U(ii)=snakposition(ii).vector(1)/80;
        V(ii)=snakposition(ii).vector(2)/80;
        str=num2str(snakposition(ii).index);
        text(X(ii)+U(ii),Y(ii)+V(ii),str)
    end
    
    
    
end

function []=PlotSnaxelLoop(figh,axh,snakposition,snaxel)
    % Plots the snaxels as arrows on the plot
    figure(figh)
    axes(axh)
    snaxInd=[snaxel(:).index];
    for jj=1:length(snaxel)
        line=[snaxel(jj).index,snaxel(jj).snaxnext];
        for ii=1:length(line)
            currSnaxSub=FindObjNum(snakposition,line(ii),snaxInd);
            X(ii)=snakposition(currSnaxSub).coord(1);
            Y(ii)=snakposition(currSnaxSub).coord(2);
        end
        plot(X,Y,'ro--')
    end

end

function []=PlotSnaxelLoopDir(figh,axh,snakposition,snaxel)
    % Plots the snaxels as arrows on the plot
    figure(figh)
    axes(axh)
    snaxInd=[snaxel(:).index];
    for jj=1:length(snaxel)
        line=[snaxel(jj).index,snaxel(jj).snaxnext];
        for ii=1:length(line)
            currSnaxSub=FindObjNum(snakposition,line(ii),snaxInd);
            X(ii)=snakposition(currSnaxSub).coord(1);
            Y(ii)=snakposition(currSnaxSub).coord(2);
            
        end
        U=X(2)-X(1);
        
        V=Y(2)-Y(1);
        quiver(X(1),Y(1),U,V,0)
    end

end

function []=PlotVolFrac(figh,axh,coord,frac)
    figure(figh)
    axes(axh)
    if frac==0
        text(coord(:,1),coord(:,2),num2str(frac),'HorizontalAlignment','center')
    else
    text(coord(:,1),coord(:,2),num2str(frac,'%.1e'),'HorizontalAlignment','center')
    end
    hold on
end

%% Forcing Function 1
% Including snaxel velocity calculation and edge normals

function [snaxel]=VelocityCalculationExpand(snaxel,snakposition)
    % simply sets the velocity of the snaxel to 1 (to be replaced with the
    % required forcing function)
    
    for ii=1:length(snaxel)
        if snaxel(ii).isfreeze
            velocity=0;
            velMultiplier=0;
        else
            velocity=1;
            norm1=snakposition(ii).normvector{1};
            norm2=snakposition(ii).normvector{2};
            dirVec=snakposition(ii).vector;
            normalDirVec=CalcNormVec2DClockWise(dirVec);
            
            isSameSide=(dot(normalDirVec,norm1)*dot(normalDirVec,norm1))>=0;
            % vector away from vertex along one of the edges
            nextSnax=snaxel(ii).connectivity(1);
            subNSnax=FindObjNum(snakposition,nextSnax);
            tanVecEdge=snakposition(subNSnax).coord-snakposition(ii).coord;
            isConcave=(dot(tanVecEdge,norm1)+dot(tanVecEdge,norm2))<=0;
            v=[];
            v(1)=1/dot(norm1,dirVec);
            v(2)=1/dot(norm2,dirVec);
            v(isinf(v))=[];
            velMultiplier=1; % default value
            if ~isempty(v)
                if isConcave
                    if isSameSide
                        velMultiplier=min(v);
                    else
                        velMultiplier=1;
                    end
                else
                    velMultiplier=max(v);
                end
                if isinf(velMultiplier)
                    warning('infinite velocity')
                end
            end
        end
        %velMultiplier=1;
        snaxel(ii).v=velocity*velMultiplier;
    end
    
end

function [snaxel]=VelocityCalculationContract(snaxel,snakposition)
    % simply sets the velocity of the snaxel to 1 (to be replaced with the
    % required forcing function)
    
    for ii=1:length(snaxel)
        if snaxel(ii).isfreeze
            velocity=0;
            velMultiplier=0;
        else
            velocity=-1;
            norm1=snakposition(ii).normvector{1};
            norm2=snakposition(ii).normvector{2};
            dirVec=snakposition(ii).vector;
            normalDirVec=CalcNormVec2DClockWise(dirVec);
            
            isSameSide=(dot(normalDirVec,norm1)*dot(normalDirVec,norm1))>=0;
            % vector away from vertex along one of the edges
            nextSnax=snaxel(ii).connectivity(1);
            subNSnax=FindObjNum(snakposition,nextSnax);
            tanVecEdge=snakposition(subNSnax).coord-snakposition(ii).coord;
            isConcave=(dot(tanVecEdge,norm1)+dot(tanVecEdge,norm2))<=0;
            v=[];
            v(1)=1/dot(norm1,dirVec);
            v(2)=1/dot(norm2,dirVec);
            v(isinf(v))=[];
            velMultiplier=1; % default value
            if ~isempty(v)
                if isConcave
                    if isSameSide
                        velMultiplier=min(v);
                    else
                        velMultiplier=1;
                    end
                else
                    velMultiplier=max(v);
                end
                if isinf(velMultiplier)
                    warning('infinite velocity')
                end
            end
        end
        %velMultiplier=1;
        snaxel(ii).v=velocity*velMultiplier;
    end
    
end

%% Contour Normal Calculation

function [snakposition]=SnaxelNormal2(snaxel,snakposition)
    % Calculates the normal at the Snaxel (According to snakes with topology control)
    

    [contourStruct]=ContourNormal2(snaxel,snakposition);
    nSnax=length(snaxel);
    
    ind1=[contourStruct(:).index1];
    ind2=[contourStruct(:).index2];
    
    for ii=1:nSnax
        contourVecSub1=FindObjNum([],snaxel(ii).index,ind1);
        contourVecSub2=FindObjNum([],snaxel(ii).index,ind2);
        contourVecSub=[contourVecSub1,contourVecSub2];
        contourVecSub(contourVecSub==0)=[];
        
        snakposition(ii).normvector=...
            {contourStruct(contourVecSub).vector};
        
        indexMat=[[contourStruct(contourVecSub).index1]',...
            [contourStruct(contourVecSub).index2]'];
        isPrecConn=logical(sum(snaxel(ii).snaxprec==indexMat,2));
        
        isNextConn=logical(sum(snaxel(ii).snaxnext==indexMat,2));
        
        snakposition(ii).vectornext=contourStruct(contourVecSub(isNextConn)).vector;
        snakposition(ii).vectorprec=contourStruct(contourVecSub(isPrecConn)).vector;

    end
    
end

function [contourStruct]=ContourNormal2(snaxel,snakposition)
    % Calculates the normal of the contour
    
    snaxIndPos=[snakposition(:).index];
    [contourStruct]=ExtractContourSnaxel(snaxel);
    % Normal Case
    for ii=1:length(contourStruct)
        indices=[contourStruct(ii).index1,contourStruct(ii).index2];
        
        [contourStruct(ii).vector]=ContourEdgeNormal(snakposition,snaxel,...
            indices,snaxIndPos);
    end
    % Backup Method
    ind1=[contourStruct(:).index1];
    ind2=[contourStruct(:).index2];
    for ii=1:length(contourStruct)
        if sum(abs(contourStruct(ii).vector))==0
            
            [contourStruct(ii).vector]=NormalContourAlternateMethods(snakposition,...
                ii,contourStruct,ind1,ind2,snaxIndPos);
        end
    end
    for ii=1:length(contourStruct)
            [contourStruct(ii).vector]=contourStruct(ii).vector/...
                norm(contourStruct(ii).vector);
    end
    
end

function [contourStruct]=ExtractContourSnaxel(snaxel)
    % Extracts the snaxel connected for each  contour edge
    snaxInd=[snaxel(:).index]';
    snaxConnects=vertcat(snaxel(:).connectivity);
    snaxContour1=[snaxInd,snaxConnects(:,1)];
    snaxContour2=[snaxInd,snaxConnects(:,2)];
    snaxContour=[snaxContour1;snaxContour2];
    cellSimilar=FindIdenticalVector(snaxContour);
    
    arraySimilar=vertcat(cellSimilar{:});
    snaxContour=snaxContour(arraySimilar(:,1),:);
    for ii=length(snaxContour(:,1)):-1:1
        contourStruct(ii).index1=snaxContour(ii,1);
        contourStruct(ii).index2=snaxContour(ii,2);
    end
end

function [normalVector]=ContourEdgeNormal(snakposition,snaxel,indices,snaxIndPos)
    % Calculates the normal of the contour
    
    snaxSubPos=FindObjNum(snaxel,indices,snaxIndPos);
    
    [normalVector]=NormalContourBaseMethods(snakposition,snaxSubPos);
    baseVectors=vertcat(snakposition(snaxSubPos).vector);
    [normalVector]=TestNormalOutPointing(baseVectors,normalVector);
    
end

function [normalVectorOut]=TestNormalOutPointing(baseVectors,normalVector)
    % Tests that the normal points out in provided a vector pointing out is
    % specified
    for jj=1:length(baseVectors(:,1))
        testDirNV(jj)=dot(baseVectors(jj,:),normalVector);
    end
    
    if sum(testDirNV<0)
        normalVectorOut=-normalVector;
    elseif sum(testDirNV==0)==2
        normalVectorOut=[0 0];
    else
        normalVectorOut=normalVector;
    end
    
end

function [normalVector]=NormalContourBaseMethods(snakposition,snaxSubPos)
    
    coord1=snakposition(snaxSubPos(1)).coord;
    coord2=snakposition(snaxSubPos(2)).coord;
    tanVec=coord1-coord2;
    % Normal Cases for calculation of 
    if sum(abs(tanVec))~=0
        % if the tangential vector is not the 0 vector
        normalVector=CalcNormVec2DClockWise(tanVec);
    else
        % else use the sum of direction vectors
        normalVector=snakposition(snaxSubPos(1)).vector+...
            snakposition(snaxSubPos(2)).vector;
    end
    
    
end

function [normalVector]=NormalContourAlternateMethods(snakposition,currentCont,...
        contourStruct,ind1,ind2,snaxIndPos)
    % Alternate method relying on adjacent edge contours to define the
    % normal vector
    
    indices=[contourStruct(currentCont).index1,contourStruct(currentCont).index2];
    otherCont1=FindObjNum([],indices,ind1);
    otherCont2=FindObjNum([],indices,ind2);
    otherCont1(otherCont1==currentCont)=[];
    otherCont2(otherCont2==currentCont)=[];
    otherCont=[otherCont1',otherCont2'];
    otherCont(otherCont==0)=[];
    adjacentVectors=vertcat(contourStruct(otherCont).vector);
    outPointVec=sum(adjacentVectors);
    snaxSubPos=FindObjNum([],indices,snaxIndPos);
    
    coord1=snakposition(snaxSubPos(1)).coord;
    coord2=snakposition(snaxSubPos(2)).coord;
    tanVec=coord1-coord2;
    % Normal Cases for calculation of 
    if sum(abs(tanVec))~=0
        % if the tangential vector is not the 0 vector
        normalVector=CalcNormVec2DClockWise(tanVec);
    else
        % else use the sum of direction vectors
        normalVector=CalcNormVec2DClockWise(snakposition(snaxSubPos(1)).vector);
    end
    [normalVector]=TestNormalOutPointing(outPointVec,normalVector);
    
end

function normalVector=CalcNormVec2DClockWise(tanVector)
    % Calculates a vector normal to another in 2D by rotating it 90 degrees
    % clockwise
    
%     normTan=norm(tanVector);
%     if tanVector(2)~=0
%         normalVector(1)=1;
%         normalVector(2)=-normalVector(1)*tanVector(1)/tanVector(2);
%         normalVector=normalVector/norm(normalVector)*normTan;
%     elseif tanVector(1)~=0
%         normalVector(2)=1;
%         normalVector(1)=-normalVector(2)*tanVector(2)/tanVector(1);
%         normalVector=normalVector/norm(normalVector)*normTan;
%     else
%         normalVector=[0 0];
%     end
    
    sizTanVec=size(tanVector);
    tanVector=reshape(tanVector,[2 1]);
    rot90=[0 1;-1 0];
    normalVector=rot90*tanVector;
    normalVector=reshape(normalVector,sizTanVec);
    
%     if sum(abs(normalVector-normalVector2))>10^-15
%         if sum(abs(normalVector+normalVector2))>10^-15
%             warning('whyyyy?')
%         end
%     end
end

%% Velocity Calculation (External Function)

function [snaxel,snakposition,snaxelmodvel]=VelocityCalculationVolumeFraction...
        (snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    
%     [snaxel,snakposition,snaxelmodvel]=...
%         VelocityForce(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
%     
%     [snaxel,snakposition,snaxelmodvel]=...
%         VelocityNormalisedForce(snaxel,snakposition,volumefraction,coeffstructure,forceparam);

    [snaxel,snakposition,snaxelmodvel]=...
        VelocityAreaOnly(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
end
%% Volume fraction calculation

function [volumefraction,coeffstruct,cellCentredGrid]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo)
    % Calculates teh volume fraction in the old cells
    
    
    insideContourInfoIndex=[refinedGrid.edge(logical(insideContourInfo)).index];
    [cellCentredGrid]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredGrid,snakposition);
    
    for ii=1:length(cellCentredGrid)
        [cellCentredGrid(ii).areaBlock,cellCentredGrid(ii).filledvolume]=...
            ExtractCellFillInformation(cellCentredGrid(ii),insideContourInfoIndex);
    end
    [volumefraction]=ExtractVolumeFractions(cellCentredGrid,volfracconnec);
    [coeffstruct,cellCentredGrid]=VolumeFractionDifferentiated(cellCentredGrid);
end

function [volumefraction]=ExtractVolumeFractions(cellCentredGrid,volfracconnec)
    % calculates the volume fraction for each cell of the old grid
    
    volumefraction=volfracconnec.cell;
    newEdgeIndices=[cellCentredGrid(:).index];
    
    for ii=1:length(volumefraction)
        
        
        newSubs=FindObjNum(cellCentredGrid,...
            volumefraction(ii).newCellInd,newEdgeIndices);
        volumefraction(ii).splitfraction=[cellCentredGrid(newSubs).filledvolume];
        volumefraction(ii).splitvolume=[cellCentredGrid(newSubs).volume];
        volumefraction(ii).totalfraction=sum(volumefraction(ii).splitfraction);
        volumefraction(ii).totalvolume=sum(volumefraction(ii).splitvolume);
        volumefraction(ii).volumefraction=volumefraction(ii).totalfraction/...
            volumefraction(ii).totalvolume;
        if volumefraction(ii).volumefraction-1>10^-12 || volumefraction(ii).volumefraction<-10^-12
            warning('volumefraction>1')
        end
        
        
    end
    
    
end

function [cellCentredGrid]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredGrid,snakposition)
    % Extracts the snaxel data and matches it to the cells
    global unstructglobal
    snakPosInd=[snakposition(:).index];
    edgeInd=[refinedGrid.edge(:).index];
    cellInd=[cellCentredGrid(:).index];
    cellCentredGrid(1).snaxel=struct([]);
    vertIndex=unstructglobal.vertex.index;
    vertCoord=unstructglobal.vertex.coord;
    
    snakPosFields={'coord','vector','vectorprec','vectornext','normvector','edgelength'};
    for ii=1:length(snaxel)
        snaxEdge=snaxel(ii).edge;
        snaxEdgeSub=FindObjNum(refinedGrid.edge,snaxEdge,edgeInd);
        snaxCells=refinedGrid.edge(snaxEdgeSub).cellindex;
        snaxCellsSub=FindObjNum(cellCentredGrid,snaxCells,cellInd);
        snaxPosSub=FindObjNum(snakposition,snaxel(ii).index,snakPosInd);
        
        for jj=1:length(snaxCellsSub)
            iToVert=vertCoord(find(vertIndex==snaxel(ii).tovertex),:); %#ok<FNDSB> % extract vertex coordinates
            iFromVert=vertCoord(find(vertIndex==snaxel(ii).fromvertex),:); %#ok<FNDSB>
            
            snaxelCell=snaxel(ii);
            for kk=1:length(snakPosFields)
                snaxelCell.(snakPosFields{kk})=0;
            end
            
            snaxelCell.coord=snakposition(snaxPosSub).coord;
            snaxelCell.vector=snakposition(snaxPosSub).vector;
            snaxelCell.vectorprec=snakposition(snaxPosSub).vectorprec;
            snaxelCell.vectornext=snakposition(snaxPosSub).vectornext;
            snaxelCell.normvector=snakposition(snaxPosSub).normvector;
            snaxelCell.edgelength=norm(iToVert-iFromVert);
            
            cellCentredGrid(snaxCellsSub(jj)).snaxel=...
                [cellCentredGrid(snaxCellsSub(jj)).snaxel,snaxelCell];
        end
    end
    
end

function [areablock,volume]=ExtractCellFillInformation(cellStruct,insideContourInfoIndex)
    % Extracts the contour information to calculate the area
    
    nBordBlocks=length(cellStruct.snaxel)/2;
    if nBordBlocks>0
        [edgeSnak]=ExtractCellSnaxelConnectedPairs(nBordBlocks,cellStruct);
        [areablock]=ExtractBorderStructure(cellStruct,edgeSnak,nBordBlocks);
        [volume,areablock]=CalculateCellVolume(areablock,cellStruct.volume);
    else
        areablock.border=[0,1];
        isFullCell=prod(FindObjNum([],[cellStruct.edge(:).index],...
            insideContourInfoIndex)>0);
        if isFullCell
            volume=cellStruct.volume;
        else
            volume=0;
        end
    end
    
end

function [bordstruct]=BorderStructure(bordLength,bordCentre,bordNormal)
    
    
    bordstruct.length=bordLength;
    bordstruct.centre=bordCentre;
    bordstruct.normal=bordNormal;
    
end

function [bordstruct,edgeSnakSub]=AddSnaxelBorders(cellStruct,edgeSnakSub)
    % Calculates snaxel-snaxel border
    % in CCW order
    
    isCCWOrder=cellStruct.snaxel(edgeSnakSub(1)).snaxnext==cellStruct.snaxel(edgeSnakSub(2)).index;
    if ~isCCWOrder
        edgeSnakSub=edgeSnakSub(2:-1:1);
    end
    
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    
    normalVector=cellStruct.snaxel(edgeSnakSub(1)).vectornext;
    normNormalVec=norm(normalVector);
    if normNormalVec~=0
        normalVector=normalVector/norm(normalVector);
    end
    bordstruct=BorderStructure(norm(activeCoord(1,:)-activeCoord(2,:)),...
        mean(activeCoord),normalVector);
    [edgeSnakSub]=CalculateCWdirectionEdge(activeCoord,...
        normalVector,edgeSnakSub);
end

function [bordStruct,activeFromVertex,activeSnaxEdge]=...
        AddEdgeSnaxelBorders(cellStruct,snakSub,vertInd,edgeInd)
    % Creates teh Border structure for Edge snaxel pairs
    
    activeSnaxEdge=cellStruct.snaxel(snakSub).edge;
    activeSnaxEdgeSub=FindObjNum(cellStruct.edge,activeSnaxEdge,edgeInd);
    activeFromVertex=cellStruct.snaxel(snakSub).fromvertex;
    activeFromVertexSub=FindObjNum(cellStruct.vertex,activeFromVertex,vertInd);
    
    actCoord(1,:)=cellStruct.snaxel(snakSub).coord;
    actCoord(2,:)=cellStruct.vertex(activeFromVertexSub).coord;
    
    bordStruct=BorderStructure(norm(actCoord(1,:)-actCoord(2,:)),...
        mean(actCoord),cellStruct.edge(activeSnaxEdgeSub).normalvector);
    
    
    
end

function [bordstruct]=AddEdgeBorders(activeStartVertex,previousEdge,cellStruct,edgeSnakSub...
        ,edgeInd,vertEdgeInd)
    % Adds the borders which are related to edges
    
    
    ll=1;
    while activeStartVertex~=cellStruct.snaxel(edgeSnakSub).fromvertex
        
        previousEdgeSub=FindObjNum(cellStruct.edge,previousEdge,edgeInd);
        
        
        nextEdgeLog=logical(sum(vertEdgeInd==activeStartVertex,2));
        nextEdgeLog(previousEdgeSub)=false;
        nextEdge=cellStruct.edge(nextEdgeLog).index;
        
        bordstruct(ll)=BorderStructure(...
            cellStruct.edge(nextEdgeLog).edgelength,...
            cellStruct.edge(nextEdgeLog).edgecentre,...
            cellStruct.edge(nextEdgeLog).normalvector);
        ll=ll+1;
        
        
        nextEdgeVertices=cellStruct.edge(nextEdgeLog).vertexindex;
        nextActiveVertex=nextEdgeVertices(nextEdgeVertices~=activeStartVertex);
        previousEdge=nextEdge;
        activeStartVertex=nextActiveVertex;
    end
    
    if ~exist('bordstruct','var');bordstruct=struct([]);end;
end

function [areablock]=ExtractBorderStructure(cellStruct,edgeSnak,nBordBlocks)
    % Extracts the border information and sets it out in the right
    % structure form
    vertInd=[cellStruct.vertex(:).index];
    edgeInd=[cellStruct.edge(:).index];
    vertEdgeInd=vertcat(cellStruct.edge(:).vertexindex);
    
    edgeSnakSub=FindObjNum(cellStruct.snaxel,edgeSnak(:,1));
    edgeSnakSub(:,2)=FindObjNum(cellStruct.snaxel,edgeSnak(:,2));
    for ii=nBordBlocks:-1:1

        % Calculate snaxel to snaxel border
        [bordstruct1,edgeSnakSub(ii,:)]=AddSnaxelBorders(cellStruct,edgeSnakSub(ii,:));
        % First snax-edge bord
        [bordstruct2,activeFromVertex,activeSnaxEdge]=...
            AddEdgeSnaxelBorders(cellStruct,edgeSnakSub(ii,1),vertInd,edgeInd);
        % Edge Borders
        activeStartVertex=activeFromVertex;
        previousEdge=activeSnaxEdge;
        bordstructEdges=AddEdgeBorders(activeStartVertex,previousEdge,...
            cellStruct,edgeSnakSub(ii,2),edgeInd,vertEdgeInd);
        
        % First snax-edge bord
        [bordstruct3,~,~]=...
            AddEdgeSnaxelBorders(cellStruct,edgeSnakSub(ii,2),vertInd,edgeInd);
        areablock(ii).blockstruct=[bordstruct1,bordstruct2,bordstructEdges,bordstruct3];
    end
    
end

function [edgeSnak]=ExtractCellSnaxelConnectedPairsProblem(nBordBlocks,cellStruct)
    % Extract the connected snaxels from cell information
    snaxConn=vertcat(cellStruct.snaxel(:).connectivity);
    snaxInd=[cellStruct.snaxel(:).index];
    snaxIndWorking=snaxInd;
    snaxConnWorking=snaxConn;
    edgeSnak=zeros([nBordBlocks,2]);
    for ii=1:nBordBlocks
        connLog=false(size(snaxIndWorking));
        for jj=1:2
            connLog=connLog | (snaxConnWorking(1,jj)==snaxIndWorking);
        end
        connLog(1)=true;
        edgeSnak(ii,:)=snaxIndWorking(connLog)';
        snaxIndWorking(connLog)=[];
        connRmvSub=find(connLog);
        snaxConnWorking(connRmvSub,:)=[]; %#ok<FNDSB>
    end
    if ~isempty(snaxIndWorking)
        warning('You''re fucking up')
    end
    
end

function [edgeSnak]=ExtractCellSnaxelConnectedPairs(nBordBlocks,cellStruct)
    % Extract the connected snaxels from cell information
    snaxNext=vertcat(cellStruct.snaxel(:).snaxnext);
    snaxInd=[cellStruct.snaxel(:).index];

    edgeSnak=zeros([nBordBlocks,2]);
    kk=0;
    for ii=1:length(snaxInd)
        nextIsInCell=sum(snaxNext(ii)==snaxInd)>0;
        if nextIsInCell
            kk=kk+1;
            edgeSnak(kk,:)=[snaxInd(ii),snaxNext(ii)];
           
        end
    end
    
    % test statements
    testEdgeSnak=numel(RemoveIdenticalEntries(edgeSnak(:)))~=numel(RemoveIdenticalEntries(snaxInd));
    if kk~=nBordBlocks || testEdgeSnak
        warning('You''re fucking up')
    end
    
end

function [edgeSnakSub]=CalculateCWdirectionEdge(activeCoord,normalVector,edgeSnakSub)
    
    testVector=activeCoord-(ones(2,1)*mean(activeCoord));
    baseVector=normalVector;
    [vecAngles]=ExtractAngle360(baseVector,testVector);
    [~,orderedInd]=sort(vecAngles);
    edgeSnakSub=edgeSnakSub(orderedInd);
end

function [volume]=CalculateGreensVolume(bordstruct)
    % THis function uses green's theorem to calculate the area inside a set
    % of snaxels
    volumeArray=zeros(size(bordstruct));
    for ii=1:length(bordstruct)
        volumeArray(ii)=dot(bordstruct(ii).centre,bordstruct(ii).normal)...
            *bordstruct(ii).length;
    end
    volume=sum(volumeArray)/2;
    
end

function [volume,areablock]=CalculateCellVolume(areablock,totalVol)
    % Calculate the volume bounded by snaxels within a cell
    
    for ii=1:length(areablock)
        areablock(ii).volume=CalculateGreensVolume(areablock(ii).blockstruct);
    end
    volume=sum(([areablock(:).volume]));
    
    if volume>totalVol
        volume=totalVol+sum(([areablock(:).volume]-totalVol));
    end
    
    
end

%% Volume Derivative calculation

function [coeffstruct,cellCentredGrid]=VolumeFractionDifferentiated(cellCentredGrid)
    % Calculates the volume coefficients for the derivative of the volume fraction
    % In order to calculate the velocities
    
    for ii=1:length(cellCentredGrid)
        [cellCentredGrid(ii).coeffBlock]=...
            ExtractCellDerivInformation(cellCentredGrid(ii));
    end
    
    coeffstruct=[cellCentredGrid(:).coeffBlock];
    for ii=1:length(coeffstruct)
        
        coeffstruct(ii).value=CalculateVolDerivCoeff(coeffstruct(ii));
    end
    
end

function [coeffblock]=ExtractCellDerivInformation(cellStruct)
    % Extracts the contour information to calculate the area
    
    nBordBlocks=length(cellStruct.snaxel)/2;
    if nBordBlocks>0
        [edgeSnak]=ExtractCellSnaxelConnectedPairs(nBordBlocks,cellStruct);
        [coeffblock]=ExtractBorderDerivStructure(cellStruct,edgeSnak,nBordBlocks);
        
    else
        coeffblock=struct([]);
    end
    
end

function [areablock]=ExtractBorderDerivStructure(cellStruct,edgeSnak,nBordBlocks)
    % Extracts the informaion for the calculation of the derivative of the
    % Area
    vertInd=[cellStruct.vertex(:).index];
    edgeInd=[cellStruct.edge(:).index];
    
    edgeSnakSub=FindObjNum(cellStruct.snaxel,edgeSnak(:,1));
    edgeSnakSub(:,2)=FindObjNum(cellStruct.snaxel,edgeSnak(:,2));
    
    for ii=nBordBlocks:-1:1
        % Calculate snaxel to snaxel border
        [bordstruct1,edgeSnakSub(ii,:)]=AddSnaxelDeriv(cellStruct,edgeSnakSub(ii,:));
        % First snax-edge border
        [bordstruct2,~,~]=...
            AddEdgeSnaxelDeriv(cellStruct,edgeSnakSub(ii,1),vertInd,edgeInd,2);
        % First snax-edge border
        [bordstruct3,~,~]=...
            AddEdgeSnaxelDeriv(cellStruct,edgeSnakSub(ii,2),vertInd,edgeInd,1);
        
        areablock(ii).blockstruct=[bordstruct1,bordstruct2,bordstruct3];
    end
    areablock=[areablock(:).blockstruct];
    
end

function [derivstruct,edgeSnakSub]=AddSnaxelDeriv(cellStruct,edgeSnakSub)
    % Calculates snaxel-snaxel border
    
    isCCWOrder=cellStruct.snaxel(edgeSnakSub(1)).snaxnext==cellStruct.snaxel(edgeSnakSub(2)).index;
    if ~isCCWOrder
        edgeSnakSub=edgeSnakSub(2:-1:1);
    end
    
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    tanVector=activeCoord(2,:)-activeCoord(1,:);
    normalVector=cellStruct.snaxel(edgeSnakSub(1)).vectornext;
    
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    normVecLength=normalVector*norm(tanVector);

    for ii=2:-1:1
        jj=ii;%abs(ii-3);
        derivstruct(ii)=...
            BorderDerivCoeffStructure(cellStruct.index,cellStruct.snaxel(edgeSnakSub(ii)).index,...
            cellStruct.snaxel(edgeSnakSub(ii)).vector*cellStruct.snaxel(edgeSnakSub(ii)).edgelength...
            ,mean(activeCoord),normVecLength,jj);
    end
    
    
end

function [derivStruct,activeFromVertex,activeSnaxEdge]=...
        AddEdgeSnaxelDeriv(cellStruct,snakSub,vertInd,edgeInd,snaxelBorderOrder)
    % Creates the Border structure for Edge snaxel pairs
    
    activeSnaxEdge=cellStruct.snaxel(snakSub).edge;
    activeSnaxEdgeSub=FindObjNum(cellStruct.edge,activeSnaxEdge,edgeInd);
    activeFromVertex=cellStruct.snaxel(snakSub).fromvertex;
    activeFromVertexSub=FindObjNum(cellStruct.vertex,activeFromVertex,vertInd);
    
    actCoord(1,:)=cellStruct.snaxel(snakSub).coord;
    actCoord(2,:)=cellStruct.vertex(activeFromVertexSub).coord;
    
    normVecLength=cellStruct.edge(activeSnaxEdgeSub).normalvector...
        *norm(actCoord(1,:)-actCoord(2,:));
    
    
    derivStruct=...
        BorderDerivCoeffStructure(cellStruct.index,cellStruct.snaxel(snakSub).index,...
        cellStruct.snaxel(snakSub).vector*cellStruct.snaxel(snakSub).edgelength,...
        mean(actCoord),normVecLength,snaxelBorderOrder);
    
end

function [derivstruct]=BorderDerivCoeffStructure(cellIndex,snaxelIndex,Dg,bordCentre,bordNormal,ordercorner)
    
    derivstruct.cellindex=cellIndex;
    derivstruct.snaxelindex=snaxelIndex;
    derivstruct.diffgrid=Dg;
    derivstruct.centre=bordCentre;
    derivstruct.normal=bordNormal;
    derivstruct.ordercorner=ordercorner; % Defines wether is negative or positive in coeff cal
    
end

function [value]=CalculateVolDerivCoeff(coeffstruct)
    % calculates the value of the velocity derivative coefficient
    
    DgVec=coeffstruct.diffgrid;
    normalVec=coeffstruct.normal;
    positionVec=coeffstruct.centre;
    orderInd= coeffstruct.ordercorner;
    rotNeg90=[0 -1; 1 0];
    
    value=(dot(DgVec,normalVec)/2)...
        +(((-1)^(orderInd-1))*dot((rotNeg90*DgVec'),positionVec));
    
end


%% Cell centred information
function [cellCentredGrid]=CellCentredGridInformation(refinedGrid)
    % Returns cell centred information about the grid being used
    
    cellCentredGrid=refinedGrid.cell;
    edgeCellInfo=vertcat(refinedGrid.edge(:).cellindex);
    
    origEdgeIndex=[refinedGrid.edge(:).index];
    origVertexIndex=[refinedGrid.vertex(:).index];
    origCellIndex=[refinedGrid.cell(:).index];
    
    for ii=1:length(cellCentredGrid)
        edgeCellLog=sum((edgeCellInfo==cellCentredGrid(ii).index),2)>0;
        totEdges=sum(edgeCellLog);
        cellCentredGrid(ii).edge(1:totEdges)=refinedGrid.edge(edgeCellLog);
        cellVertex=[cellCentredGrid(ii).edge(:).vertexindex];
        cellVertex=RemoveIdenticalEntries(cellVertex');
        cellVertexSub=FindObjNum(refinedGrid.vertex,cellVertex,origVertexIndex);
        totVertices=length(cellVertex);
        cellCentredGrid(ii).vertex(1:totVertices)=refinedGrid.vertex(cellVertexSub);
    end
    
    [cellCentredGrid]=CalculateEdgeCellNormals(cellCentredGrid);
    [cellCentredGrid]=CalculateStartCellVolumes(cellCentredGrid);
end

function [cellCentredGrid]=CalculateStartCellVolumes(cellCentredGrid)
    
    for ii=1:length(cellCentredGrid)
        
        for jj=1:length(cellCentredGrid(ii).edge)
            
            cellCentredGrid(ii).gridareablock.blockstruct(jj).centre=...
                cellCentredGrid(ii).edge(jj).edgecentre;
            cellCentredGrid(ii).gridareablock.blockstruct(jj).normal=...
                cellCentredGrid(ii).edge(jj).normalvector;
            cellCentredGrid(ii).gridareablock.blockstruct(jj).length=...
                cellCentredGrid(ii).edge(jj).edgelength;
            
        end
        [cellCentredGrid(ii).volume,cellCentredGrid(ii).gridareablock]=...
            CalculateCellVolume(cellCentredGrid(ii).gridareablock,0);
        
    end
    
    
end

function [cellCentredGrid]=CalculateEdgeCellNormals(cellCentredGrid)
    
    for ii=1:length(cellCentredGrid)
        blockEdges=vertcat(cellCentredGrid(ii).edge(:).vertexindex);
        [cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(blockEdges);
        orderedVertices=cellOrderedVertex{1}(:,1);
        orderedEdges=[cellCentredGrid(ii).edge(cellOrderedEdges{1}).index];
        orderedVertexSub=FindObjNum(cellCentredGrid(ii).vertex,orderedVertices);
        coord=vertcat(cellCentredGrid(ii).vertex(orderedVertexSub).coord);
        [isCCW]=CCWLoop(coord);
        for jj=1:length(orderedEdges)
            
            edgeVertSub=FindObjNum(cellCentredGrid(ii).vertex,cellOrderedVertex{1}(jj,:));
            edgeTanVector=cellCentredGrid(ii).vertex(edgeVertSub(2)).coord...
                -cellCentredGrid(ii).vertex(edgeVertSub(1)).coord;
            normalVector=CalcNormVec2DClockWise(edgeTanVector);
            normalVector=normalVector/norm(normalVector);
            
            [vecAngles]=ExtractAnglepm180(edgeTanVector,normalVector);
            if (isCCW && vecAngles>0) || (~isCCW && vecAngles<0)
                normalVector=-normalVector;
            end
            cellCentredGrid(ii).edge(cellOrderedEdges{1}(jj)).normalvector=normalVector;
            vertCoord=vertcat(cellCentredGrid(ii).vertex(edgeVertSub).coord);
            cellCentredGrid(ii).edge(cellOrderedEdges{1}(jj)).edgecentre=...
                mean(vertCoord);
            cellCentredGrid(ii).edge(cellOrderedEdges{1}(jj)).edgelength=...
                norm(vertCoord(1,:)-vertCoord(2,:));
        end
        
    end
    
    
end

function [volfracconnec]=VolumeFractionConnectivity(oldGrid,...
        connectGrid,cellCentreGrid,refinedGrid)
    % Modifies the connectivity information to match what is required for
    % the calculation and use of the volume fraction information
    [volfracconnec.cell]=VolumeFractionConnectivityCell(oldGrid,...
        connectGrid,cellCentreGrid);
    [volfracconnec.edge]=VolumeFractionConnectivityEdge(oldGrid,...
        volfracconnec.cell,refinedGrid);
end

function [volfracconnec]=VolumeFractionConnectivityCell(oldGrid,connectGrid,cellCentreGrid)
    % creates the right connectivity information for cells
    oldCellInd=[oldGrid.cell(:).index];
    newCellInd=[cellCentreGrid(:).index];
    splitCellInd=[connectGrid.cell(:).old];
    
    for ii=1:length(oldCellInd)
        volfracconnec(ii).oldCellInd=(oldGrid.cell(ii).index);
        splitSub=FindObjNum([],volfracconnec(ii).oldCellInd,splitCellInd);
        if splitSub==0
            newSub=FindObjNum([],volfracconnec(ii).oldCellInd,newCellInd);
            if newSub==0
                error('There is an issue with the grid splitting match')
            end
            newCellOldIndex=cellCentreGrid(newSub).index;
        else
            newCellOldIndex=connectGrid.cell(splitSub).new;
        end
        volfracconnec(ii).newCellInd=newCellOldIndex;
        volfracconnec(ii).targetfill=oldGrid.cell(ii).fill;
    end
    
end

function [volfracconnec]=VolumeFractionConnectivityEdge(oldGrid,connectCell,refinedGrid)
    % Creates the right connectivity information for edges
    oldCellInd=[oldGrid.cell(:).index];
    newEdgeInd=[refinedGrid.edge(:).index];
    subConnCellVec=[];
    for ii=1:length(connectCell)
        subConnCellVec=[subConnCellVec,ones(1,length(connectCell(ii).newCellInd))];
    end
    connectCellNewCell=[connectCell(:).newCellInd];
    
    connectCellNewCellMat=connectCellNewCell'*[1 1];
    
    for ii=length(newEdgeInd):-1:1
        
        volfracconnec(ii).newedge=newEdgeInd(ii);
        volfracconnec(ii).newcell=refinedGrid.edge(ii).cellindex;
        edgeNewCell=volfracconnec(ii).newcell;
        oldCell=zeros(size(edgeNewCell));
        oldCell(edgeNewCell==0)=0;
        
        edgeNewCellMat=ones([length(connectCellNewCell),1])*edgeNewCell;
        kk=sum(connectCellNewCellMat==edgeNewCellMat,2)>0;
%         for jj=1:length(edgeNewCell)
%             if edgeNewCell(jj)~=0
%                
%                 for kk=1:length(connectCell)
%                     if sum(connectCell(kk).newCellInd==edgeNewCell(jj))>0
%                         oldNewCell=connectCell(kk).oldCellInd;
%                     end
%                     
%                 end
%                 oldCell(jj)=oldNewCell;
%             end
%         end
        oldCell(edgeNewCell~=0)=[connectCell(subConnCellVec(kk)).oldCellInd];
        volfracconnec(ii).oldCell=oldCell;
    end
    
end


%% Snaxel position 
% Increments position

function [snakposition]=PositionSnakes(snaxel,unstructured)
    % Returns an array with Snaxel coordinates preceded by snaxel indices
    vertIndex=unstructured.vertex.index;
    vertCoord=unstructured.vertex.coord;
    fromVertex=[snaxel(:).fromvertex];
    toVertex=[snaxel(:).tovertex];
    
    nSnaxel=length(snaxel);
    
    for ii=nSnaxel:-1:1
        iToVert=vertCoord(find(vertIndex==toVertex(ii)),:); %#ok<FNDSB> % extract vertex coordinates
        iFromVert=vertCoord(find(vertIndex==fromVertex(ii)),:); %#ok<FNDSB>
        
        snakposition(ii).index=snaxel(ii).index;
        snakposition(ii).coord=iFromVert+(iToVert-iFromVert)*snaxel(ii).d;
        snakposition(ii).vectornotnorm=(iToVert-iFromVert);
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

function snaxel=SnaxelDistanceUpdate(snaxel,dt,dtSnax,maxDist)
    % Updates the current distance of the snaxel
    
    for ii=1:length(snaxel)
        movDist=snaxel(ii).v*dt;
        
        if abs(movDist)>abs(maxDist(ii))
            movDist=maxDist(ii);
        end
        snaxel(ii).d=movDist+snaxel(ii).d;
    end
    
end

function [finishedSnakesSub]=ArrivalCondition(snaxel)
    % Calculates the arrival condition for repopulation
    
    global arrivalTolerance
    
    isFreeze=[snaxel(:).isfreeze]; % finds all unfrozen
    isFwd=[snaxel(:).v]>0;
    isArrived=[snaxel(:).d]>=(1-arrivalTolerance);
    finishedSnakesSub=find((isFwd & ~isFreeze & isArrived));
end

function [insideContourInfo]=UpdateInsideContourInfo(insideContourInfo,newInsideEdges,snaxel)
    % Updates the logical inside contour info
    
    global unstructglobal
    
    edgeInd=unstructglobal.edge.index;
    snaxEdges=[snaxel(:).edge];
    [snaxEdges]=RemoveIdenticalEntries(snaxEdges);
    
    stillSnaxEdge=FindObjNum([],newInsideEdges,snaxEdges)>0;
    newInsideEdges(stillSnaxEdge)=[];
    
    newInsideEdgesSub=FindObjNum([],newInsideEdges,edgeInd);
    insideContourInfo=logical(insideContourInfo);
    insideContourInfo(newInsideEdgesSub)=true;
    
end

%% Snaxel Repopulation

function [snaxel,insideContourInfo]=...
        SnaxelRepopulate(unstructured,snaxel,insideContourInfo)
    
    global unstructglobal
    % adds snaxels once a corner has been reached
    [finishedSnakes]=ArrivalCondition(snaxel);
    edgeSnaxel=[snaxel(:).edge];
    edgeVertIndex=unstructured.edge.vertexindex;
    edgeIndex=unstructured.edge.index;
    
    %newInsideEdges=[snaxel(finishedSnakes).edge];
    %insideContourInfo(newInsideEdges)=1;
    
    [snaxel,newInsideEdges,delIndex]=RepopIterativeBreedingProcess...
        (snaxel,finishedSnakes,edgeVertIndex,edgeIndex,edgeSnaxel,unstructglobal);
    
    % Removing from repopulation list snaxels that would hit a vertex
    % which has already bred
    [delIndex]=RemoveIdenticalEntries(delIndex(:));
    snaxel=DeleteSnaxel(snaxel,delIndex);
    
    
    [insideContourInfo]=UpdateInsideContourInfo(insideContourInfo,...
        newInsideEdges,snaxel);
    
    
end

function [snaxel,newInsideEdges,delIndex]=RepopIterativeBreedingProcess...
        (snaxel,finishedSnakes,edgeVertIndex,edgeIndex,edgeSnaxel,unstructglobal)
    % Iterative Breeding process for the breeding of edges
    
    newInsideEdges=[];
    kk=0;
    delIndex=[];
    while ~isempty(finishedSnakes)
        kk=kk+1;
        breedSub=finishedSnakes(1);
        finishedSnakes(1)=[];
        newInsideEdges(kk)=snaxel(breedSub).edge; %#ok<AGROW>
        %insideContourInfo(newInsideEdgesSub(kk))=1;
        snaxelIndexStart=max([snaxel(:).index])+1;
        % Doubling the snaxel to be removed
        [snaxel,indexDoubled,snaxelIndexStart]=...
            RepopAddDuplicateSnaxel(snaxel,snaxelIndexStart,breedSub);
        % Find CCW order of two same snax
        [connec]=RepopExtractConnectionOrder(snaxel,breedSub,indexDoubled,unstructglobal);
        
        % Generate the new Snaxels
        snaxelRepop=InitialSnaxelStructure(snaxel(breedSub).tovertex,edgeVertIndex,...
            edgeIndex,snaxel(breedSub).edge,snaxelIndexStart,...
            connec);
        [snaxel]=AddSnaxel(snaxel,snaxelRepop);
        % Remove breeding snaxels snaxels
        delIndex(kk,1:2)=[snaxel(breedSub).index,indexDoubled]; %#ok<AGROW>
        
        for jj=1:length(snaxelRepop)
            rmvFinSnakes=(edgeSnaxel(finishedSnakes)==snaxelRepop(jj).edge);
            finishedSnakes(rmvFinSnakes)=[];
        end
    end
    
end

function [snaxel,indexDoubled,snaxelIndexStart]=...
        RepopAddDuplicateSnaxel(snaxel,snaxelIndexStart,ii)
    additionalsnaxel=snaxel(ii);
    additionalsnaxel.index=snaxelIndexStart;
    connecRemove=additionalsnaxel.connectivity(1);
    additionalsnaxel.connectivity(1)=snaxel(ii).index;
    isNextRmv=additionalsnaxel.snaxnext==connecRemove;
    if isNextRmv
        additionalsnaxel.snaxnext=snaxel(ii).index;
    else
        additionalsnaxel.snaxprec=snaxel(ii).index;
    end
    
    indexDoubled=snaxelIndexStart;
    snaxelIndexStart=snaxelIndexStart+1;
    [snaxel]=AddSnaxel(snaxel,additionalsnaxel);
end

function [connec]=RepopExtractConnectionOrder(snaxel,ii,indexDoubled,unstructglobal)
    % orders the connectivity indices in CCW order for the repopulation stage
    snaxelIndices=[snaxel(:).index];
    [posSnaxSub]=FindObjNum(snaxel,[snaxel(ii).index,indexDoubled],snaxelIndices);
    [posSnaxSub1]=FindObjNum(snaxel,[snaxel(posSnaxSub).connectivity],snaxelIndices);
    [posSnaxSub2]=FindObjNum(snaxel,[snaxel(posSnaxSub1).connectivity],snaxelIndices);
    [posSnaxSub3]=FindObjNum(snaxel,[snaxel(posSnaxSub2).connectivity],snaxelIndices);
    posSnaxSub=[posSnaxSub',posSnaxSub1',posSnaxSub2',posSnaxSub3'];
    posSnaxSub=RemoveIdenticalEntries(posSnaxSub);
    [snakposition]=PositionSnakes(snaxel(posSnaxSub),unstructglobal);
    
    [precSnax,nextSnax]=CCWNeighbours(snaxel,snaxel(ii).index,snakposition);
    if precSnax==indexDoubled
        connec=[indexDoubled,snaxel(ii).index];
    elseif nextSnax==indexDoubled
        connec=[snaxel(ii).index,indexDoubled];
    end
end

function [snaxel]=AddSnaxel(snaxel,additionsnaxel)
    % [snaxel]=AddSnaxel(snaxel,additionsnaxel)
    % adds a snaxel defined in the standard way and connects it up in between
    % the specified connection snaxels
    
    [connection,indexSnaxCon]=ExtractConnection(additionsnaxel);
    %connection=additionsnaxel.connectivity;
    subConnection=FindObjNum(snaxel,connection);
    for ii=1:length(connection)
        % in the iith connected snaxel find the index to the other connection
        singlesnaxel=snaxel(subConnection(ii));
        connecReplace=indexSnaxCon(ii);
        connecRemove=connection(connection~=connection(ii));
        [singlesnaxel]=ModifyConnection(singlesnaxel,connecRemove,connecReplace);
        snaxel(subConnection(ii))=singlesnaxel;
    end
    snaxel=[snaxel,additionsnaxel];
end

function [singlesnaxel]=ModifyConnection(singlesnaxel,connecRemove,connecReplace)
        % Removes a connection to a snaxel from another snaxel and replaces
        % it with the replacement connection
        refresh=singlesnaxel.connectivity==connecRemove;
        isNext=singlesnaxel.snaxnext==connecRemove;
        isPrec=singlesnaxel.snaxprec==connecRemove;
        if sum(refresh)==0
            error('invalid connection to break');
        end
        if isNext && isPrec
            warning('Topology Collapsing');
        end
        % then replace that number to the snaxel that is being introduced
        singlesnaxel.connectivity(refresh)=connecReplace;
        if isNext
            singlesnaxel.snaxnext=connecReplace;
        elseif isPrec
            singlesnaxel.snaxprec=connecReplace;
        end
end


%% Time Update

function maxDist=MaxTravelDistance(snaxel)
    % Calculates the maximum distance that can be travelled by a snaxel
    dSnax=[snaxel(:).d];
    vSnax=[snaxel(:).v];
    edgeSnax=[snaxel(:).edge];
    nSnax=length(edgeSnax);
    
    maxDist=ones(1,nSnax);
    
    for ii=1:nSnax
        sameEdgeSnax=find(edgeSnax==edgeSnax(ii));
        sameEdgeSnax(sameEdgeSnax==ii)=[];
        
        if vSnax(ii)>=0
            maxDist(ii)=1-dSnax(ii);
        else
            maxDist(ii)=-dSnax(ii);
        end
        if numel(sameEdgeSnax)==1
            impactDist=(1-dSnax(ii)-dSnax(sameEdgeSnax))/...
                (vSnax(ii)+vSnax(sameEdgeSnax))*vSnax(ii);
            if vSnax(ii)+vSnax(sameEdgeSnax)==0
                impactDist=1;
            end
            
            if vSnax(ii)>=0 && impactDist>=0
                maxDist(ii)=min([maxDist(ii),impactDist]);
            elseif vSnax(ii)<0 && impactDist<=0
                maxDist(ii)=max([maxDist(ii),impactDist]);
            end
            
        elseif numel(sameEdgeSnax)>1
            warning('More than 2 snaxels on the same edge')
        end
    end
    
end

function [dt,dtSnax2,maxDist]=TimeStepCalculation(snaxel,maxStep,maxDt,dtMin)
    % Calculates the timestep to ensure that no snaxel crosses the line without
    % crossing an intersection
    
    
    canMove=[snaxel(:).isfreeze]==0;
    maxDist=MaxTravelDistance(snaxel);
    vSnax=[snaxel(:).v];
    dtSnax=maxDist./vSnax;
    dtSnax2=dtSnax;
    dtSnax2(~canMove)=0;
    dtSnax2(dtSnax2<0)=0;
    dtSnax2(dtSnax2>maxDt)=maxDt;
    dtSnax2(isnan(dtSnax2))=0;
    
    dtSnax=dtSnax(canMove);
    dtSnax=dtSnax(dtSnax>=0);
    dtStep=min(abs(maxStep./vSnax(canMove)));
    dtMax=min(dtSnax);
    dt=min([dtStep, dtMax, maxDt]);
    if dt==0
        warning('Time Step was forced to minimum set by user')
        dt=dtMin;
    end
    if dt<0
        warning('Time going back')
    end
end

%% Order Snaxels
% Functions which order snaxels in a CCW manner
function connecOrder=CCWOrderAroundNode(snaxelRepop,baseEdge)
    % Orders snaxel leaving from a single node in counter clockwise order
    
    global unstructglobal
    
    [snaxRepopPos]=PositionSnakes(snaxelRepop,unstructglobal);
    connecSnax=[snaxelRepop(:).index];
    nSnaxRepop=length(snaxelRepop);
    connecSub=1:length(snaxelRepop);
    
    baseEdgeSub=FindObjNum([],baseEdge,unstructglobal.edge.index);
    baseVertexIndex=unstructglobal.edge.vertexindex(baseEdgeSub,:);
    baseVertexSub=FindObjNum([],baseVertexIndex,unstructglobal.vertex.index);
    
    fromVertLog=(snaxelRepop(1).fromvertex==baseVertexIndex);
    toVertLog=(snaxelRepop(1).fromvertex~=baseVertexIndex);
    
    fromCoord=unstructglobal.vertex.coord(baseVertexSub(fromVertLog),:);
    toCoord=unstructglobal.vertex.coord(baseVertexSub(toVertLog),:);
    baseVector=toCoord-fromCoord;
    
    
    testCoord=zeros(nSnaxRepop,2);
    testVector=zeros(nSnaxRepop,2);
    
    for ii=1:length(connecSub)
        testCoord(ii,:)=snaxRepopPos(connecSub(ii)).coord;
        testVector(ii,:)=snaxRepopPos(connecSub(ii)).vector;
    end
    
    %[vecAngles]=ExtractAnglepm180(baseVector,testVector);
    [vecAngles]=ExtractAngle360(baseVector,testVector);
    [~,chainOrder]=sort(vecAngles);
    
    connecOrder=connecSnax(chainOrder);
end

function [connection,indexSnaxCon]=ExtractConnection(additionsnaxel)
    % finds the snaxels which are included in the connectivity information but
    % are not part of the set of snaxel presented
    
    nSnax=length(additionsnaxel);
    additionIndex=[additionsnaxel(:).index];
    connectivity=zeros(nSnax,2);
    for ii=1:nSnax
        connectivity(ii,:)=additionsnaxel(ii).connectivity;
    end
    for ii=1:nSnax
        rmvCon=connectivity==additionIndex(ii);
        if sum(rmvCon(:))~=0
            connectivity(rmvCon)=0;
        end
        
    end
    [subSnaxCon,yCol]=find(connectivity);
    indexSnaxCon=additionIndex(subSnaxCon); % returns the index of the edge of the snaxel block
    for ii=1:length(subSnaxCon)
        connection(ii)=connectivity(subSnaxCon(ii),yCol(ii)); % returns the snaxel outside the block
    end
    
    
    
end

function [precSnax,nextSnax]=CCWNeighbours(snaxel,snakInd,snakposition)
    % gives the counter clockwise order of neighbouring points
    snaxelIndices=[snaxel(:).index];
    snaxelIndicesPos=[snakposition(:).index];
    snakSub=FindObjNum(snaxel,snakInd,snaxelIndices);
    snakSubP=FindObjNum(snakposition,snakInd,snaxelIndicesPos);
    
    
    precSnax=snaxel(snakSub).snaxprec;
    nextSnax=snaxel(snakSub).snaxnext;
    
end

%% Cleaning conquest

function [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo)
    % Function containing the process allowing the deletion of Snaxels at each
    % time step
    delIndex=0;
    kk=0;
    
    while ~isempty(delIndex)
        canMove=[snaxel(:).isfreeze]==0;
        if ~isempty(snaxel(canMove))
            [delIndex]=FindBadSnaxels(snaxel(canMove),insideContourInfo);
        else
            delIndex=[];
        end
        delSub=FindObjNum(snaxel,delIndex);
        insideEdgesInd=[snaxel(delSub).edge];
        snaxel=DeleteSnaxel(snaxel,delIndex);
        [insideContourInfo]=UpdateInsideContourInfo(insideContourInfo,...
            insideEdgesInd,snaxel);
        kk=kk+1;
    end
    
    
end

function [delIndex]=FindBadSnaxels(snaxel,insideContourInfo)
    % Function containing all the deleting conditions for the Snaxels
    delIndex=[];
    %[wanderingIndex]=FindWanderingSnaxels(snaxel);
     insideIndex=[];
    %[insideIndex]=FindInsideSnaxels(snaxel,insideContourInfo);
   
    [edgeConIndex]=FindEdgeConnectedSnaxels(snaxel);
    
    %delIndex=[insideIndex;edgeConIndex;arrivedIndex];
    delIndex=[insideIndex;edgeConIndex];
    [delIndex]=RemoveIdenticalEntries(delIndex);
end

function [delIndex]=FindInsideSnaxels(snakes,insideContourInfo)
    % Finds Snaxels on an edge inside the contour
    delIndex=[];
    global unstructglobal
    edgeInd=unstructglobal.edge.index;
    edgeSnakes=[snakes(:).edge];
    
    edgeSub=FindObjNum([],edgeSnakes,edgeInd);
    
    for ii=1:length(snakes)
        
        if insideContourInfo(edgeSub(ii))
            delIndex=[delIndex;snakes(ii).index];
        end
    end
    
end

function [delIndex]=FindEdgeConnectedSnaxels(snaxel)
    % Finds Snaxels on an edge inside the contour
    
    nSnakes=length(snaxel);
    
    snaxelEdgeVec=[snaxel(:).edge];
    for ii=1:nSnakes
        snaxelConnec1Vec(ii)=snaxel(ii).connectivity(1);
        snaxelConnec2Vec(ii)=snaxel(ii).connectivity(2);
    end
    snaxelIndex=[snaxel(:).index];
    
    [edgeMatVert,edgeMatHoriz]=meshgrid(snaxelEdgeVec,snaxelEdgeVec);
    snaxelConnec1Mat=ones(nSnakes,1)*snaxelConnec1Vec;
    snaxelConnec2Mat=ones(nSnakes,1)*snaxelConnec2Vec;
    snaxelIndexMat=snaxelIndex'*ones(1,nSnakes);
    
    equalEdge=edgeMatVert==edgeMatHoriz; % Finds snakes on the same edge
    connectEdge=(snaxelIndexMat==snaxelConnec2Mat)|...
        (snaxelIndexMat==snaxelConnec1Mat); % Finds snakes that are connected
    problemEdge=equalEdge & connectEdge & (~eye(nSnakes)); % Finds snakes which are both
    
    xProb= sum(problemEdge); % remove the snakes which are both
    delIndex=[snaxel(logical(xProb)).index]';
    
    
end

function snaxel=DeleteSnaxel(snaxel,delIndex)
    
    numDel=length(delIndex);
    snaxelIndices=[snaxel(:).index];
    delSnaxSub=zeros([1,numDel]);
    for ii=1:numDel
        delSnaxSub(ii)=FindObjNum(snaxel,delIndex(ii),snaxelIndices);
        
        connect=snaxel(delSnaxSub(ii)).connectivity;
        targSnaxSub=FindObjNum(snaxel,connect,snaxelIndices);
        
        
        for jj=1:length(targSnaxSub)
            singlesnaxel=snaxel(targSnaxSub(jj));
            connecReplace=connect(connect~=connect(jj));
            connecRemove=delIndex(ii);
            if ~isempty(connecReplace)
                [singlesnaxel]=ModifyConnection(singlesnaxel,connecRemove,connecReplace);
                snaxel(targSnaxSub(jj))=singlesnaxel;
            else
                warning('Check Topology Collapse')
            end
            
        end
        
    end
    
    snaxel(delSnaxSub)=[];

end

function [insideContourInfoRev]=ReverseInsideContourInfo(snaxel,insideContourInfo,unstructured)
    % Revereses the information in insideContourInfo
    
    edgeInd=unstructured.edge.index;
    snaxEdgeSub=FindObjNum([],[snaxel(:).edge],edgeInd);
    
    insideContourInfoRev=~insideContourInfo;
    insideContourInfoRev(snaxEdgeSub)=false;
    
    
    
    
end


%% Snaxel freezing functions

function [snaxel]=FreezingFunction(snaxel,borderVertices,mergeTopo)
    % Function which selects snaxels that need to be frozen
    
    if ~exist('mergeTopo','var'),mergeTopo=false;end
    
    if mergeTopo
        edgeFreeze=[];
    else
        [edgeFreeze]=FreezeEdgeContact(snaxel);
    end
    
    [borderFreeze]=FreezeBorderContact(snaxel,borderVertices);
    freezeIndex=[edgeFreeze,borderFreeze];
    
    [freezeIndex]=RemoveIdenticalEntries(freezeIndex);
    if ~isempty(freezeIndex)
        freezeSub=FindObjNum(snaxel,freezeIndex);
        for ii=1:length(freezeIndex)
            snaxel(freezeSub(ii)).isfreeze=1;
        end
    end
end

function [freezeIndex,pairs]=FreezeEdgeContact(snaxel,saveSingleVal)
    % Function which returns the index of Snaxels which have met in the middle
    % of an edge
    
    if ~exist('saveSingleVal','var');saveSingleVal=false;  end
    pairs=[];
    
    dSnax=[snaxel(:).d];
    vSnax=[snaxel(:).v];
    isFreeze=[snaxel(:).isfreeze];
    indexSnax=[snaxel(:).index];
    edgeSnax=[snaxel(:).edge];
    fromvertSnax=[snaxel(:).fromvertex];
    nSnax=length(edgeSnax);
    
    isImpact=false(1,nSnax);
    
    for ii=1:nSnax
        sameEdgeSnax=find(edgeSnax==edgeSnax(ii));
        sameEdgeSnax(sameEdgeSnax==ii)=[];
        
        
        if numel(sameEdgeSnax)>1
            warning('More than 2 snaxels on the same edge')
            
            sameEdgeDeltaD=dSnax(sameEdgeSnax)-dSnax(ii);
            sameEdgeDeltaDBelow=sameEdgeDeltaD;
            sameEdgeDeltaDAbove=sameEdgeDeltaD;
            sameEdgeDeltaDBelow(sameEdgeDeltaDBelow>0)=-2;
            sameEdgeDeltaDAbove(sameEdgeDeltaDAbove<0)=2;
            [~,prevSnax]=max(sameEdgeDeltaDBelow);
            [~,nextSnax]=min(sameEdgeDeltaDAbove);
            
            
        elseif numel(sameEdgeSnax)==1
            isImpact(ii)=EdgeImpactCondition(dSnax,vSnax,fromvertSnax,ii,sameEdgeSnax);
            if saveSingleVal && isImpact(ii) && ~isFreeze(ii)
                pairs=[indexSnax(ii),indexSnax(sameEdgeSnax)];
                break
            end
            if isFreeze(ii)
                isImpact(ii)=false;
            end
            
        end
    end
    freezeIndex=indexSnax(isImpact);
end

function [isImpact]=EdgeImpactCondition(dSnax,vSnax,fromvertSnax,sub1,sub2)
    % edge Impact condition
    
    global arrivalTolerance
    sameDir=fromvertSnax(sub1)==fromvertSnax(sub2);
    if sameDir
        warning('Snaxels are adjacent and moving in same direction Impact condition is invalid')
        dSnax(sub2)=1-dSnax(sub2);
        vSnax(sub2)=-vSnax(sub2);
    end
    
    deltaD=(1-dSnax(sub1)-dSnax(sub2));
    if deltaD>=0
        signDeltaD=1;
    else
        signDeltaD=-1;
    end
    isContact=deltaD==0;
    deltaV=signDeltaD*(vSnax(sub1)+vSnax(sub2));
    
    isClose=abs(deltaD)<=arrivalTolerance;
    isApproaching=(deltaV)>0;
    isImpact=(isClose && isApproaching) || isContact;
end

function [freezeIndex]=FreezeBorderContact(snaxel,borderVertices)
    % Function which returns the index of Snaxels which have hit the edge of
    % the design space
    
    
    [arrivedSub]=ArrivalCondition(snaxel);
    arrivedIndex=[snaxel(arrivedSub).index];
    
    indexSnax=[snaxel(arrivedSub).index];
    destVertex=[snaxel(arrivedSub).tovertex];
    nSnax=length(destVertex);
    
    isImpact=false(1,nSnax);
    
    for ii=1:nSnax
        borderDestination=sum(borderVertices==destVertex(ii));
        
        if borderDestination
            isImpact(ii)=true;
        end
        
    end
    freezeIndex=indexSnax(isImpact);
    
end

function [borderVertices]=FindBorderVertex(unstructured)
    
    edgeVertIndex=unstructured.edge.vertexindex;
    edgeIndex=unstructured.edge.index;
    vertexIndex=unstructured.vertex.index;
    edgeCellIndex=unstructured.edge.cellindex;
    
    nVert=length(vertexIndex);
    
    isBorderVertex=false([nVert 1]);
    ii=1;
    vertex=0;
    while ~isempty(vertex)
        vertex(1)=[];
        edges=FindEdgesIndex(vertexIndex(ii),edgeVertIndex,edgeIndex);
        vertEdgeSub=FindObjNum([],edges,edgeIndex);
        vertCellIndex=edgeCellIndex(vertEdgeSub,:);
        isOnBorderCell=sum(sum((vertCellIndex==0)))>0;
        if isOnBorderCell && ~isBorderVertex(ii)
            isBorderVertex(ii)=true;
            vertex=edgeVertIndex(vertEdgeSub,:);
            vertex(vertex(:)==vertexIndex(ii))=[];
        end
        
        if ~isempty(vertex)
            ii=find(vertex(1)==vertexIndex);
        end
    end
    
    borderVertices=vertexIndex(isBorderVertex);
end

function [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo)
    % function which merges topologies
    
    global unstructglobal
    edgeInd=unstructglobal.edge.index;
    
    workingPair=[];
    [mergeIndex,workingPair]=FreezeEdgeContact(snaxel,true);
    ll=1;
    savMergeInd=[];
    snaxOrig=snaxel;
    while ~isempty(workingPair)
        
        
        for ii=1:length(workingPair)
            [precSnax(ii),nextSnax(ii)]=CCWNeighbours(snaxel,workingPair(ii),snakposition);
        end
        if ~isempty(workingPair)
            [snaxel]=MergeTopologies(snaxel,workingPair,precSnax,nextSnax);
            [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
            savMergeInd(ll)=mergeIndex;
            ll=ll+1;
        end
        clear workingPair precSnax nextSnax
        [mergeIndex,workingPair]=FreezeEdgeContact(snaxel,true);
    end
    snaxInd=[snaxOrig(:).index];
    insideSnaxSub=FindObjNum(snaxel,savMergeInd,snaxInd);
    insideEdgesInd=[snaxOrig(insideSnaxSub).edge];
    
    [insideContourInfo]=UpdateInsideContourInfo(insideContourInfo,...
        insideEdgesInd,snaxel);
    
end

function [snaxel]=MergeTopologies(snaxel,workPair,precSnak,nextSnak)
    % Merges topologies by inverting connections
    
    workSub=FindObjNum(snaxel,workPair);
    snakSave=snaxel(workSub);
    snaxel=DeleteSnaxel(snaxel,workPair);
    
    snaxInd=[snaxel(:).index];
    precSub=FindObjNum(snaxel,precSnak,snaxInd);
    nextSub=FindObjNum(snaxel,nextSnak,snaxInd);
    
    
    for ii=1:length(workPair)
        jj=abs(ii-3);
        snaxel(precSub(ii)).snaxnext=nextSnak(jj);
        snaxel(precSub(ii)).connectivity=[snaxel(precSub(ii)).snaxnext,...
            snaxel(precSub(ii)).snaxprec];
        snaxel(nextSub(ii)).snaxprec=precSnak(jj);
        snaxel(nextSub(ii)).connectivity=[snaxel(nextSub(ii)).snaxnext,...
            snaxel(nextSub(ii)).snaxprec];
        
    end
    
end

function [snaxel]=MergeTopologies2(snaxel,workPair,precSnak,nextSnak)
    % Merges topologies by inverting connections
    
    workSub=FindObjNum(snaxel,workPair);
    snakSave=snaxel(workSub);
    snaxel=DeleteSnaxel(snaxel,workPair);
    
    snaxInd=[snaxel(:).index];
    precSub=FindObjNum(snaxel,precSnak,snaxInd);
    nextSub=FindObjNum(snaxel,nextSnak,snaxInd);
    
    
    for ii=1:length(workPair)
        jj=abs(ii-3);
        
        connecPrec=snaxel(precSub(ii)).connectivity;
        connecNext=snaxel(nextSub(jj)).connectivity;
        
        connecPrec(connecPrec==nextSnak(ii))=nextSnak(jj);
        connecNext(connecNext==precSnak(jj))=precSnak(ii);
        
        snaxel(precSub(ii)).connectivity=connecPrec;
        snaxel(nextSub(jj)).connectivity=connecNext;
    end
    
end

%% Various

function [loopsnaxel]=OrderSurfaceSnaxel(snaxel)
    % function extracting the snaxels into their separate loops
    global unstructglobal
    
    snaxPositions=PositionSnakes(snaxel,unstructglobal);
    
    nSnax=length(snaxel);
    blockSegments=zeros(2*nSnax,2);
    for ii=1:nSnax
        for jj=0:1
            blockSegments(2*ii-jj,:)=[snaxel(ii).index,snaxel(ii).connectivity(jj+1)];
        end
    end
    
    cellSimilar=FindIdenticalVector(blockSegments);
    for ii=1:length(cellSimilar)
        blockEdgeIndex(ii)=cellSimilar{ii}(1);
    end
    blockEdges=blockSegments(blockEdgeIndex,:);
    % Order edges into closed loops
    [cellOrderedVertex]=OrderBlockEdges(blockEdges);
    snaxIndex=[snaxel(:).index];
    for ii=1:length(cellOrderedVertex)
        loopsnaxel(ii).snaxel.index=[cellOrderedVertex{ii}(:,1);cellOrderedVertex{ii}(1:2,1)];
        loopIndices=FindObjNum(snaxel,loopsnaxel(ii).snaxel.index,snaxIndex);
        loopsnaxel(ii).snaxel.coord=vertcat(snaxPositions(loopIndices).coord);
        %loopsnaxel(ii).edge.index=isEdgeIndex(cellOrderedEdges{ii});
    end
    
end

function [snaxelrev]=ReverseSnakes(snaxel)
    % Reverses a set of snaxels such that they are pointing the other way
    
    for ii=length(snaxel):-1:1
        snaxelrev(ii)=snaxel(ii);
        % unchanged
        snaxelrev(ii).index=snaxel(ii).index;
        snaxelrev(ii).edge=snaxel(ii).edge;
        snaxelrev(ii).isfreeze=snaxel(ii).isfreeze;
        snaxelrev(ii).connectivity=snaxel(ii).connectivity;
        % changed
        snaxelrev(ii).d=1-snaxel(ii).d;
        snaxelrev(ii).v=-snaxel(ii).v;
        snaxelrev(ii).tovertex=snaxel(ii).fromvertex;
        snaxelrev(ii).fromvertex=snaxel(ii).tovertex;
        snaxelrev(ii).snaxprec=snaxel(ii).snaxnext;
        snaxelrev(ii).snaxnext=snaxel(ii).snaxprec;
    end
    
    
    
end

function [snaxelrev]=ReverseSnakesConnection(snaxel)
    % Reverses a set of snaxels such that they are pointing the other way
    
    for ii=length(snaxel):-1:1
        snaxelrev(ii)=snaxel(ii);
        % changed
        snaxelrev(ii).snaxprec=snaxel(ii).snaxnext;
        snaxelrev(ii).snaxnext=snaxel(ii).snaxprec;
    end
    
    
    
end

%% Copied/Modified from main

function []=template()
    
end

%%%%%%%%%%%%%%%%%%%%%%%
%% OLD Working  CODE %%
%%%%%%%%%%%%%%%%%%%%%%%

%% Do not delete

%{
function [precSnax,nextSnax]=CCWNeighbours2(snaxel,snakInd,snakposition)
    % gives the counter clockwise order of neighbouring points
    snaxelIndices=[snaxel(:).index];
    snaxelIndicesPos=[snakposition(:).index];
    snakSub=FindObjNum(snaxel,snakInd,snaxelIndices);
    snakSubP=FindObjNum(snakposition,snakInd,snaxelIndicesPos);
    connecInd=snaxel(snakSub).connectivity;
    %connecSub=FindObjNum(snakposition,connecInd);
    baseCoord=snakposition(snakSubP).coord;
    
    [~,connecNonSimInd,intermidiateCon]=FindFirstNonSimilarSnaxel(snaxel,connecInd,snakInd,...
        baseCoord,snakposition);
    
    [contourStruct]=ExtractContourStruct(snaxel,snakposition,connecNonSimInd,...
        intermidiateCon,snaxelIndices);
    baseVec=sum(vertcat(contourStruct(:).vector));
    [precSnaxNonSim,nextSnaxNonSim]=CCWConnections(snakposition,snakInd,...
        connecNonSimInd,baseVec);
    
    precSnaxLog=precSnaxNonSim==connecNonSimInd;
    nextSnaxLog=nextSnaxNonSim==connecNonSimInd;
    
    precSnax=connecInd(precSnaxLog);
    nextSnax=connecInd(nextSnaxLog);
    
end

function [contourStruct]=ExtractContourStruct(snaxel,snakposition,connecNonSimInd,...
        intermidiateCon,snaxelIndices)
    
    for ii=1:length(connecNonSimInd)
        truncInd=[connecNonSimInd(ii),intermidiateCon{ii}(end)];
        truncSub=FindObjNum(snaxel,truncInd,snaxelIndices);
        snaxelTrunc=snaxel(truncSub);
        for jj=1:length(snaxelTrunc)
            
            snaxelTrunc(jj).connectivity=snaxelTrunc(jj).connectivity...
                (snaxelTrunc(jj).connectivity==truncInd(1) |...
                snaxelTrunc(jj).connectivity==truncInd(2))*[1 1];
            
        end
        [contourStruct(ii)]=ContourNormal2(snaxelTrunc,snakposition);
    end
end

function [precSnax,nextSnax]=CCWConnections(snakposition,snakInd,connecIndex,outBaseVec)
    % gives the counter clockwise order of neighbouring points
    snakSub=FindObjNum(snakposition,snakInd);
    connecSub=FindObjNum(snakposition,connecIndex);
    
    baseCoord=snakposition(snakSub).coord;
    %baseVec=snakposition(snakSub).vector;
    connecCoord=vertcat(snakposition(connecSub).coord);
    
    connecVec=zeros(size(connecCoord));
    for ii=1:length(connecIndex)
        connecVec(ii,:)=connecCoord(ii,:)-baseCoord;
    end
    
    vecAngles=ExtractAngle360(outBaseVec,connecVec);
    is0Angle=find(vecAngles==0, 1);
    if ~isempty(is0Angle)
        error('0 Angles should Be a Thing of the past')
        %{
        if numel(is0Angle)==1
            not0Angle=find(vecAngles~=0);
            if vecAngles(not0Angle)>pi
                vecAngles(is0Angle)=0;
            elseif vecAngles(not0Angle)<pi
                vecAngles(is0Angle)=2*pi;
            elseif vecAngles(not0Angle)==pi
                warning('Calculation of CCW connections is problematic')
            end
        else
            warning('Calculation of CCW connections is problematic')
        end
        %}
    end
    
    [~,iNext]=min(vecAngles);
    [~,iPrec]=max(vecAngles);
    precSnax=connecIndex(iPrec);
    nextSnax=connecIndex(iNext);
    
end

function [snakposition]=SnaxelNormal(snaxel,snakposition)
    % Calculates the normal at the Snaxel (According to snakes with topology control)
    
    [normalcontourvec]=ContourNormal(snaxel,snakposition);
    nSnax=length(snaxel);
    nNCV=length(normalcontourvec);
    normalVecVertex=zeros(nNCV,2);
    for ii=1:nNCV
        normalVecVertex(ii,1)=normalcontourvec(ii).vertex(1);
        normalVecVertex(ii,2)=normalcontourvec(ii).vertex(2);
    end
    
    for ii=1:nSnax
        contourVecSub=(normalVecVertex(:,1)==snaxel(ii).index);
        snakposition(ii).normvector=...
            {normalcontourvec(contourVecSub).vector};
        
        
%         if numel(find(normalVecVertex==snaxel(ii).index))>2
%             warning('normal contour vec not operating as expected')
%         end
    end
    
end

function [normalcontourvec]=ContourNormal(snaxel,snakposition)
    % Calculates the normal of the contour
    
    snaxInd=[snaxel(:).index];
    nSnax=length(snaxel);
    ll=0;
    
    
    for ii=1:nSnax
        for kk=1:2
            ll=ll+1;
            snakIndex=snakposition(ii).index;
            
            baseCoord=snakposition(ii).coord;
            
            nextSnakIndex=snaxel(ii).connectivity(kk);
            nextSnakSub=FindObjNum(snaxel,nextSnakIndex,snaxInd);
            normalcontourvec(ll).vertex=[snakIndex,nextSnakIndex];
            
            tanVec=baseCoord-snakposition(nextSnakSub).coord;
            if sum(abs(tanVec))~=0
                % if the tangential vector is not the 0 vector
                normalVector=CalcNormVec2D(tanVec);
            else
                % else use the direction vectors
                normalVector=sum([snakposition(ii).vector;...
                    snakposition(nextSnakSub).vector]);
            end
            if sum(abs(normalVector))==0
                [nextSub,nextSnax]=...
                    FindFirstNonSimilarSnaxel(snaxel,nextSnakIndex,snakIndex,...
                    baseCoord,snakposition);
                nextSnakSub=FindObjNum(snaxel,nextSnax,snaxInd);
                tanVec=baseCoord-snakposition(nextSnakSub).coord;
                if sum(abs(tanVec))~=0
                    % if the tangential vector is not the 0 vector
                    normalVector=CalcNormVec2D(tanVec);
                else
                    % else use the direction vectors
                    normalVector=sum([snakposition(ii).vector;...
                        snakposition(nextSnakSub).vector]);
                end
            end
            for jj=1:2
                workIndex=normalcontourvec(ll).vertex(jj);
                workIndex=FindObjNum(snaxel,workIndex,snaxInd);
                testDirNV(jj)=dot(snakposition(workIndex).vector,normalVector);
            end
            if sum(testDirNV<0)
                normalVector=-normalVector;
                for jj=1:2
                    workIndex=normalcontourvec(ll).vertex(jj);
                    workIndex=FindObjNum(snaxel,workIndex,snaxInd);
                    testDirNV(jj)=dot(snakposition(workIndex).vector,normalVector);
                end
                if sum(testDirNV<0)
                    %                    warning('Snaxel moving into the volume')
                end
            end
            if sum(abs(normalVector))==0
                error('vertex Normal is [0 0]')
            end
            normalcontourvec(ll).vector=normalVector/norm(normalVector);
        end
    end
end

function [connecSub,connecSnax,intermediateCon]=...
        FindFirstNonSimilarSnaxel(snaxel,connecSnax,snaxIndex,baseCoord,snaxPositions)
    % Checks the case where snaxels have the same coordinates and orders them
    % based on the following snaxel
    
    snaxInd=[snaxel(:).index];
    snaxPInd=[snaxPositions(:).index];
    connecSub=FindObjNum(snaxel,connecSnax,snaxInd);
    connecSubP=FindObjNum(snaxPositions,connecSnax,snaxPInd);
    intermediateCon{length(connecSub)}=[];
    for ii=1:length(connecSub)
        precConnectionIndex=snaxIndex; 
        ll=1;
        intermediateCon{ii}(ll)=precConnectionIndex;
        ll=ll+1;
        while sum(baseCoord~=snaxPositions(connecSubP(ii)).coord)==0;
            replaceConnec=snaxel(connecSub(ii)).connectivity;
            replaceConnec=replaceConnec(replaceConnec~=precConnectionIndex);
            precConnectionIndex=connecSnax(ii);
            intermediateCon{ii}(ll)=precConnectionIndex;
            ll=ll+1;
            connecSnax(ii)=replaceConnec;
            connecSub(ii)=FindObjNum(snaxel,connecSnax(ii),snaxInd);
            connecSubP(ii)=FindObjNum(snaxPositions,connecSnax(ii),snaxPInd);
        end
    end
    
end


function [derivstruct,edgeSnakSub]=AddSnaxelDeriv2(cellStruct,edgeSnakSub)
    % Calculates snaxel-snaxel border
    
    snaxNormVector=vertcat(cellStruct.snaxel(edgeSnakSub).normvector);
    for jj=1:2
        for ll=1:2
            if sum(snaxNormVector{1,jj}==snaxNormVector{2,ll})==2
                snaxNormVecIndex=jj;
            end
        end
    end
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    
    normVecLength=snaxNormVector{1,snaxNormVecIndex}*norm(activeCoord(1,:)-activeCoord(2,:));
    
    [edgeSnakSub]=CalculateCWdirectionEdge(activeCoord,...
        snaxNormVector{1,snaxNormVecIndex},edgeSnakSub);
    for ii=2:-1:1
        jj=ii;%abs(ii-3);
        derivstruct(ii)=...
            BorderDerivCoeffStructure(cellStruct.index,cellStruct.snaxel(edgeSnakSub(ii)).index,...
            cellStruct.snaxel(edgeSnakSub(ii)).vector*cellStruct.snaxel(edgeSnakSub(ii)).edgelength...
            ,mean(activeCoord),normVecLength,jj);
    end
    
    
end


function [bordstruct,edgeSnakSub]=AddSnaxelBorders2(cellStruct,edgeSnakSub)
    % Calculates snaxel-snaxel border
    
    snaxNormVector=vertcat(cellStruct.snaxel(edgeSnakSub).normvector);
    for jj=1:2
        for ll=1:2
            if sum(snaxNormVector{1,jj}==snaxNormVector{2,ll})==2
                snaxNormVecIndex=jj;
            end
        end
    end
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    
    bordstruct=BorderStructure(norm(activeCoord(1,:)-activeCoord(2,:)),...
        mean(activeCoord),snaxNormVector{1,snaxNormVecIndex});
    [edgeSnakSub]=CalculateCWdirectionEdge(activeCoord,...
        snaxNormVector{1,snaxNormVecIndex},edgeSnakSub);
end

%}
