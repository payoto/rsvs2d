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

% %#codegen

%% Main execution functions
function [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=SnakesSensitivity(refinedGrid,looprestart,...
        oldGrid,connectionInfo,param)
    % JUST DECLARING VARIABLES FOR LATER
    
    % plotInterval,numSteps,makeMov,boundstr
    global arrivalTolerance unstructglobal maxStep maxDt
    
    
    varExtract={'arrivalTolerance','maxStep','maxDt','snakesConsole'};
    [arrivalTolerance,maxStep,maxDt,snakesConsole]=ExtractVariables(varExtract,param);
    
    % ACTUALLY DOING STUFF
    refinedGriduns=ModifReshape(refinedGrid);
    oldGridUns=ModifReshape(oldGrid);
    
    unstructglobal=refinedGriduns;
    %profile on
    [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=...
        RunSnakesProcess(refinedGriduns,refinedGrid,looprestart,...
        oldGrid,oldGridUns,connectionInfo,param);
    %profile viewer
    
    if snakesConsole
        figure,semilogy(1:length(snakSave),[snakSave(:).currentConvVolume])
        title('Volume error')
        ylabel('Root Mean squared error on volume convergence')
        xlabel('number of iterations')
    end
end

function [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=...
        RunSnakesProcess(refinedGriduns,refinedGrid,looprestart,...
        oldGrid,oldGridUns,connectionInfo,param)
    % Main execution container for Snakes
    
    % Unpacking NECESSARY variables
    global maxStep maxDt snaxInitPos
    
    varExtract={'mergeTopo','boundstr','snakesConsole','dtRatio','snaxInitPos','checkSensitivities'};
    [mergeTopo,boundstr,snakesConsole,dtRatio,snaxInitPos,checkSensitivities]=ExtractVariables(varExtract,param);
    
    forceparam=param.snakes.force;
    
    
    % Starting process
    disp(['  Start initialisation'])
    tStepStart=now;
    [cellCentredGrid,volfracconnec,borderVertices,snaxel,...
        insideContourInfo]=InitialisationRestart(refinedGriduns,...
        refinedGrid,looprestart,oldGrid,connectionInfo,mergeTopo,boundstr,param);
   
    
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snakposition]=SnaxelNormal2(snaxel,snakposition);
    %CheckResultsLight(refinedGriduns,snakposition,snaxel)
    
    
    [freezeVertex,borderEdgInd]=IdentifyProfileEdgeVertices(refinedGrid);
    
    tStepEnd=now;
    disp(['  Initialisation time:',datestr(tStepEnd-tStepStart,'HH:MM:SS:FFF')]);
    tStart=now;
   
    tEnd=now;
    disp(['  ',int2str(ii),' Steps Performed'])
    disp(['  Iteration Time:',datestr(tEnd-tStart,'HH:MM:SS:FFF')]);
    disp(['  Volume converged to ',num2str(currentConvVolume,'%.5e')])
    
    [snaxel,snakposition,loopsnaxel]=FinishSnakes(snaxel,...
        borderVertices,refinedGriduns,param);
    
        GetSnaxelSensitivities(snaxel,refinedGriduns,refinedGrid,volfracconnec,...
                cellCentredGrid,insideContourInfo,forceparam);
    
    
    restartsnake.snaxel=snaxel;
    restartsnake.insideContourInfo=insideContourInfo;
    restartsnake.cellCentredGrid=cellCentredGrid;
    restartsnake.volfracconnec=volfracconnec;
    restartsnake.borderVertices=borderVertices;
    
    
    
end

function [cellCentredGrid,volfracconnec,borderVertices,snaxel,...
        insideContourInfo]=InitialisationRestart(refinedGriduns,...
        refinedGrid,looprestart,oldGrid,connectionInfo,mergeTopo,boundstr,param)
    
    varExtract={'restart'};
    [restart]=ExtractVariables(varExtract,param);
    if ~restart
        [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
            StartSnakeProcess(refinedGriduns,refinedGrid,looprestart,...
            oldGrid,connectionInfo,mergeTopo,boundstr);
    else
        [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
            RestartSnakeProcess(looprestart);
    end
    
    
end

function [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        StartSnakeProcess(refinedGriduns,refinedGrid,loop,...
        oldGrid,connectionInfo,mergeTopo,boundstr)
    disp('    Generate Cell Centred Grid')
    [cellCentredGrid]=CellCentredGridInformation(refinedGrid);
    disp('    Establish Cell Colume Information')
    [volfracconnec]=VolumeFractionConnectivity(oldGrid,...
        connectionInfo,cellCentredGrid,refinedGrid);
    
    
    insideContourInfo=refinedGriduns.edge.(boundstr{2});
    disp('    Find Border Vertices')
    [borderVertices]=FindBorderVertex(refinedGriduns);
    disp('    Initialise Snaxel Grid')
    
    % Inside contour info will change depending on the type of contour under
    % consideration (0 contour or 1 contour)
    [snaxel,insideContourInfo]=SnaxelInitialisation(refinedGriduns,loop,insideContourInfo,boundstr);
    
    [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo);
    [snaxel]=FreezingFunction(snaxel,borderVertices,mergeTopo);
end


function [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        RestartSnakeProcess(restartsnake)
    
    [cellCentredGrid]=restartsnake.cellCentredGrid;
    [volfracconnec]=restartsnake.volfracconnec;
    [borderVertices]=restartsnake.borderVertices;
    
    snaxel=restartsnake.snaxel;
    insideContourInfo=restartsnake.insideContourInfo;
    
end

function [snaxel,snakposition,loopsnaxel]=FinishSnakes(snaxel,...
        borderVertices,refinedGriduns,param)
    
    varExtract={'edgeFinish','TEShrink','LEShrink','axisRatio'};
    [edgeFinish,TEShrink,LEShrink,axisRatio]=ExtractVariables(varExtract,param);
    
    %disp('Finished Iterations , starting Post Process')
    [snaxel]=FreezingFunction(snaxel,borderVertices);
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    
    %disp('Creating Snaxel Loops')
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel);
    for ii=1:length(loopsnaxel)
        loopsnaxel(ii).snaxel.coord(:,2)=loopsnaxel(ii).snaxel.coord(:,2)*axisRatio;
    end
    
    switch edgeFinish
        case 'shrink'
            for ii=1:length(loopsnaxel)
                loopsnaxel(ii).snaxel.coord=ShrinkEdges(loopsnaxel(ii).snaxel.coord,LEShrink,TEShrink);
            end
        case 'sharpen'
            for ii=1:length(loopsnaxel)
                [loopsnaxel(ii).snaxel.coord]=SharpenEdges(loopsnaxel(ii).snaxel.coord,TEShrink,LEShrink);
            end
        case 'none'
            
    end
    
end


function []=GetSnaxelSensitivities(snaxel,refinedGriduns,refinedGrid,...
        volfracconnec,cellCentredGrid,insideContourInfo,forceparam)
    
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snakposition]=SnaxelNormal2(snaxel,snakposition);
    [volumefraction,coeffstructure,cellCentredGridSnax]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo);
    forceparam.isLast=true;
    forceparam.lengthEpsilon=1e-8;
    [snaxel,snakposition,snaxelmodvel,velcalcinfo,sensSnax]...
        =VelocityCalculationVolumeFraction(snaxel,snakposition,volumefraction,...
        coeffstructure,forceparam);
    testSensitivity(snaxel,snakposition,sensSnax)
end


%% Return Data

function [snakSave]=WriteSnakSave(param,snaxel,dt,snakposition,...
        volumefraction,cellCentredGridSnax,currentConvVelocity,...
        currentConvVolume,movFrame,velcalcinfo,insideContourInfo)
    
    varExtract={'snakData'};
    [snakData]=ExtractVariables(varExtract,param);
    
    switch snakData
        case 'all'
            [snakSave]=WriteSnakSaveHeavy(snaxel,dt,snakposition,...
                volumefraction,cellCentredGridSnax,currentConvVelocity,...
                currentConvVolume,movFrame,velcalcinfo,insideContourInfo);
            
        case 'light'
            [snakSave]=WriteSnakSaveLight(dt,volumefraction,currentConvVelocity,...
                currentConvVolume);
            
    end
    
    
    
    
    function [snakSave]=WriteSnakSaveHeavy(snaxel,dt,snakposition,...
            volumefraction,cellCentredGridSnax,currentConvVelocity,...
            currentConvVolume,movFrame,velcalcinfo,insideContourInfo)
        
        snakSave.snaxel=snaxel;
        snakSave.dt=dt;
        snakSave.snakposition=snakposition;
        snakSave.volumefraction=volumefraction;
        snakSave.cellCentredGrid=cellCentredGridSnax;
        snakSave.currentConvVelocity=currentConvVelocity;
        snakSave.currentConvVolume=currentConvVolume;
        snakSave.movFrame=movFrame;
        snakSave.velcalcinfo=velcalcinfo;
        snakSave.insideContourInfo=insideContourInfo;
        
    end
    
    function [snakSave]=WriteSnakSaveLight(dt,volumefraction,currentConvVelocity,...
            currentConvVolume)
        
        snakSave.dt=dt;
        
        snakSave.volumefraction.targetfill=[volumefraction(:).targetfill];
        snakSave.volumefraction.currentfraction=[volumefraction(:).volumefraction];
        snakSave.volumefraction.totVolume=[volumefraction(:).totalvolume];
        
        snakSave.currentConvVelocity=currentConvVelocity;
        snakSave.currentConvVolume=currentConvVolume;
        
        
    end
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

%% Snaxel Initialisation

function [snaxel,insideContourInfo]=SnaxelInitialisation(unstructured,loop,insideContourInfo,boundstr)
    allLoopEdgeIndex=[];
    for ii=1:length(loop)
        allLoopEdgeIndex=[allLoopEdgeIndex,loop(ii).edge.index];
    end
    for ii=1:length(loop)
        if ii==1
            snaxelIndexStart=0;
        else
            snaxelIndexStart=max([snaxel(:).index]);
        end
        [loopsnakes(ii).snaxels,insideContourInfo]=SnaxelLoop(unstructured,loop(ii),...
            snaxelIndexStart,insideContourInfo,boundstr,allLoopEdgeIndex);
        if ii==1
            snaxel=loopsnakes(ii).snaxels;
        else
            snaxel=[snaxel,loopsnakes(ii).snaxels];
        end
    end
    
end

function [snaxel,insideContourInfo]=SnaxelLoop(unstructured,loop,...
        snaxelIndexStart,insideContourInfo,boundstr,allLoopEdgeIndex)
    
    loopVertInd=loop.vertex.index(1:end-2); %RemoveIdenticalEntries(loop.vertex.index);
    loopVertSub=FindObjNum([],loopVertInd,unstructured.vertex.index);
    vertCoord=loop.vertex.coord(1:end-2,:);
    isLoopCCW=CCWLoop(vertCoord);
    if isLoopCCW
        initVertexIndex=loopVertInd;
    else
        %initVertexIndex=loop.vertex.index(end-2:-1:1);
        initVertexIndex=loopVertInd(end:-1:1);
        vertCoord=vertCoord(end:-1:1);
    end
    %initVertexIndex=loop.vertex.index(1:end-2);
    %initVertexIndex=loop.vertex.index([2:3,5:end-3]);
    edgeVertIndex=unstructured.edge.vertexindex;
    edgeIndex=unstructured.edge.index;
    loopEdgeIndex=loop.edge.index;
    vertCoordFull=unstructured.vertex.coord;
    vertIndex=unstructured.vertex.index;
    
    switch boundstr{3}
        case '1bound'
            isInside=false;
        case '0bound'
            isInside=true;
        otherwise
            error('Invalid boundstr flag')
    end
    
    
    isInside=CheckInsideFill(vertCoordFull,edgeVertIndex,initVertexIndex,...
        edgeIndex,vertIndex,insideContourInfo);
    
    [snaxel,~]=InitialSnaxelStructure(initVertexIndex,edgeVertIndex,...
        edgeIndex,loopEdgeIndex,allLoopEdgeIndex,snaxelIndexStart,isInside);
    
    [delIndex]=FindInsideSnaxels(snaxel,insideContourInfo);
    if ~isInside
        loopEdgeSubs=FindObjNum([],loopEdgeIndex,edgeIndex);
    else
        snaxInd=[snaxel(:).index];
        keepIndex=delIndex;
        keepSub=FindObjNum([],keepIndex,snaxInd);
        logDelInd=true(size(snaxInd));
        logDelInd(keepSub)=false;
        delIndex=snaxInd(logDelInd);
        loopEdgeSubs=[];
        snaxel=ReverseSnakes(snaxel);
    end
    
    snaxel=DeleteSnaxel(snaxel,delIndex);
    
    [snaxel]=TestSnaxelLoopDirection(snaxel);
    
    insideContourInfo(loopEdgeSubs)=1;
    
end

function isInside=CheckInsideFill(vertCoord,edgeVertIndex,initVertexIndex,...
        edgeIndex,vertIndex,insideContourInfo)
    
    flag=true;
    while flag
        rootVert=initVertexIndex(1);
        vertexEdges=FindEdgesIndex(rootVert,edgeVertIndex,edgeIndex);
        vertexEdgesSub=FindObjNum([],vertexEdges,edgeIndex);
        test=sum(insideContourInfo(vertexEdgesSub));
        flag=false;
        if test==0 || test==length(vertexEdgesSub)
            flag=true;
            initVertexIndex=initVertexIndex([2:end,1]);
            
            
        end
    end
    
    
    neighboursInd=sort(edgeVertIndex(vertexEdgesSub,:),2);
    
    rootEdge=[initVertexIndex(end) rootVert];
    exitEdge=[rootVert initVertexIndex(2)];
    
    [posExit]=FindMatchingVectors(neighboursInd,sort(exitEdge));
    [posRoot]=FindMatchingVectors(neighboursInd,sort(rootEdge));
    
    for ii=1:length(vertexEdges)
        neighbourVert(ii)=neighboursInd(ii,find(neighboursInd(ii,:)~=rootVert));
    end
    
    neighbourVertSub=FindObjNum([],neighbourVert,vertIndex);
    rootVertSub=FindObjNum([],rootVert,vertIndex);
    
    rootCoord=vertCoord(rootVertSub,:);
    neighCoord=vertCoord(neighbourVertSub,:);
    
    vectors=neighCoord-(ones([length(vertexEdges),1])*rootCoord);
    angles=ExtractAngle360(vectors(posRoot,:),vectors);
    
    indNoEx=1:length(angles);
    indNoEx([posExit,posRoot])=[];
    anglesExt=angles(indNoEx);
    
    indOut=indNoEx(find((anglesExt>=angles(posRoot)) & (anglesExt<=angles(posExit))));
    testIsOut=insideContourInfo(vertexEdgesSub(indOut));
    
    if sum(testIsOut)==0
        isInside=true;
    elseif sum(testIsOut)==numel(testIsOut)
        isInside=false;
    else
        error('Cannot decide if is inside or not')
    end
    
end

function [pos]=FindMatchingVectors(vecList,vec)
    
    pos=find(sum(abs(vecList-(ones([length(vecList(:,1)),1])*vec)),2)==0);
    
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
        edgeIndex,currLoopEdgeIndex,invalidEdgeIndex,snaxelIndexStart,isInside,connectivity)
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
            if numel(currLoopEdgeIndex)>1
                baseEdgeExploit=zeros(length(currLoopEdgeIndex),1);
                for jj=1:length(currLoopEdgeIndex)
                    if sum(vertexEdges==currLoopEdgeIndex(jj))
                        baseEdgeExploit(jj)=1;
                    end
                end
                baseEdgeExploitTest=baseEdgeExploit+[baseEdgeExploit(2:end);baseEdgeExploit(1)];
                baseEdgeExploit=find(baseEdgeExploitTest==2);
            else
                baseEdgeExploit=1;
            end
            
            connecOrder=CCWOrderAroundNode(snaxelNotOrdered,currLoopEdgeIndex(baseEdgeExploit));
            generateOrder=FindObjNum(snaxelNotOrdered,connecOrder);
            if isInside
                generateOrder=generateOrder(end:-1:1);
            end
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
    
    global snaxInitPos
    
    numSE=length(snaxelEdges); % provides information about snaxel from same vertex
    snaxelEdgesSub=FindObjNum([],snaxelEdges,edgeIndex);
    kkLocal=0;
    for jj=1:numSE
        kk=kk+1;
        kkLocal=kkLocal+1;
        snaxIndex=kk+snaxelIndexStart;
        cellSimVertex(jj)=snaxIndex;
        dist=snaxInitPos; % Snaxel initialisation, it starts at the vertex
        velocity=0; % Initialisation velocity
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
    snaxel.acc=0;
    snaxel.fromvertex=vertexOrig;
    snaxel.tovertex=vertexDest;
    snaxel.edge=edge;
    snaxel.connectivity=[snaxPrec snaxNext];
    snaxel.snaxprec=snaxPrec;
    snaxel.snaxnext=snaxNext;
    snaxel.isfreeze=0;
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
            vertcat(contourStruct(contourVecSub).vector);
        
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

function [snaxel,snakposition,snaxelmodvel,velcalcinfo,sensSnax]=VelocityCalculationVolumeFraction...
        (snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    velcalcinfo=[];
    sensSnax=[];
    switch forceparam.velType
        case 'default'
            [snaxel,snakposition,snaxelmodvel,velcalcinfo,sensSnax]=...
                VelocityLengthMinimisationSQP(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'force'
            [snaxel,snakposition,snaxelmodvel]=...
                VelocityForce(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'forceMin'
            [snaxel,snakposition,snaxelmodvel]=...
                VelocityForceMinimisation(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'normForce'
            [snaxel,snakposition,snaxelmodvel]=...
                VelocityNormalisedForce(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'velMinSQP'
            [snaxel,snakposition,snaxelmodvel,velcalcinfo]=...
                VelocityLengthMinimisationSQP(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'velMinLin'
            [snaxel,snakposition,snaxelmodvel,velcalcinfo]=...
                VelocityLengthMinimisation(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'velAreaOnly'
            [snaxel,snakposition,snaxelmodvel]=...
                VelocityAreaOnly(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'velForceFE'
            [snaxel,snakposition,snaxelmodvel]=...
                VelocityForce_FE(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
        case 'expand'
            [snaxel]=VelocityCalculationExpand(snaxel,snakposition);
            snaxelmodvel=snaxel;
        case 'contract'
            [snaxel]=VelocityCalculationContract(snaxel,snakposition);
            snaxelmodvel=snaxel;
            
        otherwise %'default'
            [snaxel,snakposition,snaxelmodvel,velcalcinfo,sensSnax]=...
                VelocityLengthMinimisationSQP(snaxel,snakposition,volumefraction,coeffstructure,forceparam);
    end
    
end

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
        snakposition(ii).vertInit=iFromVert;
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

function [finishedSnakesSub]=ArrivalCondition(snaxel)
    % Calculates the arrival condition for repopulation
    
    global arrivalTolerance maxDt
    v=[snaxel(:).v];
    d=[snaxel(:).d];
    isFreeze=[snaxel(:).isfreeze]; % finds all unfrozen
    isFwd=v>0;
    isArrived=d>=(1-arrivalTolerance);
    willImpact=abs((1-d)./v)<maxDt;
    finishedSnakesSub=find((isFwd & ~isFreeze & isArrived & willImpact));
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


%% Time Update

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
        canMove=[snaxel(:).isfreeze]~=1;
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

function [delIndex]=FindInsideSnaxels(snaxel,insideContourInfo)
    % Finds Snaxels on an edge inside the contour
    delIndex=[];
    global unstructglobal
    edgeInd=unstructglobal.edge.index;
    edgeSnakes=[snaxel(:).edge];
    
    edgeSub=FindObjNum([],edgeSnakes,edgeInd);
    
    for ii=1:length(snaxel)
        
        if insideContourInfo(edgeSub(ii))
            delIndex=[delIndex;snaxel(ii).index];
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

function [freezeVertex,borderEdgInd]=IdentifyProfileEdgeVertices(refinedGrid)
    
    borderEdges=[refinedGrid.edge(:).boundaryis0] & [refinedGrid.edge(:).boundaryis1];
    borderEdgInd=[refinedGrid.edge(borderEdges).index];
    borderVertices=sort([refinedGrid.edge(borderEdges).vertexindex]);
    freezeVertex=RemoveIdenticalEntries(borderVertices);
    
end


%% Check sensitivity
function []=testSensitivity(snaxel,snakposition,sensSnax)
    
    global unstructglobal
    sensSnax(:,find(sum(abs(sensSnax))==0))=[];
    dCurr=[snaxel(:).d];
    kk=1;
    snaxInd=[snaxel(:).index];
    snaxOrd=zeros(size(snaxel));
    ordList=1:numel(snaxOrd);
    for ii=1:length(snaxel)
        kk=FindObjNum([],[snaxel(kk).snaxnext],snaxInd);
        if ordList(kk)==0
            kk=min(ordList(ordList~=0));
        end
        snaxOrd(ii)=kk;
        ordList(kk)=0;
    end
    %snaxOrd(end+1)=snaxOrd(1);
    coord1=vertcat(snakposition(:).coord);
    [dir]=sum((vertcat(snakposition(:).vector)~=0).*[ones([length(snaxel), 1]),...
        ones([length(snaxel), 1])*2],2);
    [testPos1]=CreateComparisonMatrix(coord1);
    l=max(sum(vertcat(snakposition(:).vectornotnorm).^2,2));
    for ii=1:length(sensSnax(1,:)),
        snaxCopy=snaxel;
        e1=(1)./sensSnax(:,ii);
        e2=(-1)./sensSnax(:,ii);
        e1_sel=min(e1(e1>0));
        e2_sel=min(e2(e2>0));
        e_sel(ii)=min([e1_sel,e2_sel]);
        dChange{ii}=sensSnax(:,ii)/max(abs(sensSnax(:,ii)));
        %dChange{ii}=-sensSnax(:,ii)/100;
        dAct=dCurr'+dChange{ii};
        
        for jj=1:length(snaxel)
            snaxCopy(jj).d=dAct(jj);
        end
        [snakposition2]=PositionSnakes(snaxCopy,unstructglobal);
        
        
        coord2=vertcat(snakposition2(:).coord);
        [testPos2]=CreateComparisonMatrix(coord2);
        [newOrd]=CompareTestpos(testPos1,testPos2,snaxOrd,dir);
        figure
        plot(coord1(snaxOrd,1),coord1(snaxOrd,2),'o-',coord2(newOrd,1),coord2(newOrd,2),'o-')
        hold on
        for jj=1:length(newOrd)
            plot([coord1(newOrd(jj),1),coord2(newOrd(jj),1)],[coord1(newOrd(jj),2),coord2(newOrd(jj),2)],'k--')
        end
        title(['mode ',int2str(ii)])
        axis equal
        
    end
    
    
end

function []=ExecuteSensitivity(snaxel,snakposition,sensSnax)
    
    global unstructglobal
    sensSnax(:,find(sum(abs(sensSnax))==0))=[];
    dCurr=[snaxel(:).d];
    
    [snaxOrd]=SplitSnaxLoops(snaxel); % Isolates individual loops
    maxDistRatio=1/4;
    [dChange]=FindModalDistanceChange(sensSnax,maxDistRatio);
    
    
    coord1=vertcat(snakposition(:).coord);
    [dir]=sum((vertcat(snakposition(:).vector)~=0).*[ones([length(snaxel), 1]),...
        ones([length(snaxel), 1])*2],2);
    [testPos1]=CreateComparisonMatrix(coord1);
    
    l=max(sum(vertcat(snakposition(:).vectornotnorm).^2,2));
    for ii=1:length(sensSnax(1,:)),
        snaxCopy(ii,:)=snaxel;
        
        dAct=dCurr'+dChange{ii};
        for jj=1:length(snaxel)
            snaxCopy(ii,jj).d=dAct(jj);
        end
        [snakposition2(ii,:)]=PositionSnakes(snaxCopy(ii,:),unstructglobal);
    end 
        
%         coord2=vertcat(snakposition2(:).coord);
%         [testPos2]=CreateComparisonMatrix(coord2);
%         [newOrd]=CompareTestpos(testPos1,testPos2,snaxOrd,dir);
%         
        
%         figure
%         plot(coord1(snaxOrd,1),coord1(snaxOrd,2),'o-',coord2(newOrd,1),coord2(newOrd,2),'o-')
%         hold on
%         for jj=1:length(newOrd)
%             plot([coord1(newOrd(jj),1),coord2(newOrd(jj),1)],[coord1(newOrd(jj),2),coord2(newOrd(jj),2)],'k--')
%         end
%         title(['mode ',int2str(ii)])
%         axis equal
        
    
    
    
end

function [snaxOrd]=SplitSnaxLoops(snaxel)
    % Splits a snake into its component loops
    
    kk=1;
    jj=1;
    snaxInd=[snaxel(:).index];
    snaxOrd=zeros(size(snaxel));
    ordList=1:numel(snaxOrd);
    for ii=1:length(snaxel)
        kk=FindObjNum([],[snaxel(kk).snaxnext],snaxInd);
        if ordList{jj}(kk)==0
            jj=jj+1;
            kk=min(ordList(ordList~=0));
        end
        snaxOrd{jj}(ii)=kk;
        ordList(kk)=0;
    end
    
end

function [dChange]=FindModalDistanceChange(sensSnax,maxDistRatio)
    
    nMode=size(snSnax,2);
    dChange{nMode}=[];
    for ii=1:nMode
        dChange{ii}=sensSnax(:,ii)/max(abs(sensSnax(:,ii)))*maxDistRatio;
    end
end

function [testPos]=CreateComparisonMatrix(coord)
    
    [m,n]=size(coord);
    testPos{n}=[];
    for jj=1:n
        testPos{jj}=zeros(m);
        for ii=1:m
            testPos{jj}(:,ii)=(coord(ii,jj)>coord(:,jj));
        end
    end
    
    
end

function [newOrd]=CompareTestpos(t1,t2,ord,dir)
    
    for ii=1:length(t1)
        tD{ii}=t1{ii}~=t2{ii};
    end
    tX=true(size(tD{1}));
    tDel=zeros(size(tD{1}));
    for ii=1:length(tD)
        tX=tX & tD{ii};
        
    end
    for ii=1:length(tD)
        tDel=tDel+ii*(tD{ii} & ~tX);
    end
    
    newOrd=ord;
    for ii=1:length(ord)-1
        if tX(newOrd(ii+1),newOrd(ii))
            tX(newOrd(ii+1),newOrd(ii))=false;
            tX(newOrd(ii),newOrd(ii+1))=false;
            interim=newOrd(ii);
            newOrd(ii)=newOrd(ii+1);
            if ii==1
                newOrd(end)=newOrd(ii+1);
            end
            newOrd(ii+1)=interim;
        end
    end
    newOrd=newOrd(end:-1:1);
    for ii=1:length(ord)-1
        if tX(newOrd(ii+1),newOrd(ii))
            interim=newOrd(ii);
            newOrd(ii)=newOrd(ii+1);
            if ii==1
                newOrd(end)=newOrd(ii+1);
            end
            newOrd(ii+1)=interim;
        end
    end
    kk=1;
    newOrd=[newOrd(end-1),newOrd];
    rmSnak=[];
    for ii=2:length(newOrd)-1
        tTest=tDel([newOrd(ii-1),newOrd(ii),newOrd(ii+1)],[newOrd(ii)]);
        if (sum(tTest~=dir(newOrd(ii)) & tTest~=0))
            rmSnak(kk)=ii;
            kk=kk+1;
        end
    end
    newOrd(rmSnak)=[];
    
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
        loopsnaxel(ii).snaxel.index=[cellOrderedVertex{ii}(:,1)];
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
        snaxelrev(ii).acc=-snaxel(ii).acc;
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

function [points]=ShrinkEdges(points,LEShrink,TEShrink)
    % Function which allows to make sharp trailing edges and leading edges
    [~,iLE]=min(points(:,1));
    [~,iTE]=max(points(:,1));
    
    [m,~]=size(points);
    
    iLEm1=iLE-1;
    [iLEm1]=IndexMod(iLEm1,m);
    iLEp1=iLE+1;
    [iLEp1]=IndexMod(iLEp1,m);
    
    iTEm1=iTE-1;
    [iTEm1]=IndexMod(iTEm1,m);
    iTEp1=iTE+1;
    [iTEp1]=IndexMod(iTEp1,m);
    
    points(iLEm1,2)=points(iLEm1,2)-LEShrink;
    points(iLEp1,2)=points(iLEp1,2)+LEShrink;
    
    points(iTEm1,2)=points(iTEm1,2)+TEShrink;
    points(iTEp1,2)=points(iTEp1,2)-TEShrink;
    
    if points(iLEm1,2)<=points(iLE+1,2)
         points([iLEm1,iLE+1],:)=[];
    end
    if points(iTEm1,2)>=points(iTEp1,2)
         points([iTEm1,iTEp1],:)=[];
    end
    
end

function [points]=SharpenEdges(points,TEShrink,LEShrink)
    % Function which allows to make sharp trailing edges and leading edges
    [~,iLE]=min(points(:,1));
    [~,iTE]=max(points(:,1));
    
    [m,~]=size(points);
    
    iLEm1=iLE-1;
    [iLEm1]=IndexMod(iLEm1,m);
    iLEp1=iLE+1;
    [iLEp1]=IndexMod(iLEp1,m);
    iLEm2=iLE-2;
    [iLEm2]=IndexMod(iLEm2,m);
    iLEp2=iLE+2;
    [iLEp2]=IndexMod(iLEp2,m);
    
    iTEm1=iTE-1;
    [iTEm1]=IndexMod(iTEm1,m);
    iTEp1=iTE+1;
    [iTEp1]=IndexMod(iTEp1,m);
    iTEm2=iTE-2;
    [iTEm2]=IndexMod(iTEm2,m);
    iTEp2=iTE+2;
    [iTEp2]=IndexMod(iTEp2,m);
    
    if LEShrink
        [points(iLEm1,:)]=Align3points(points([iLE,iLEm2],:),points(iLEm1,:));
        [points(iLEp1,:)]=Align3points(points([iLE,iLEp2],:),points(iLEp1,:));
    end
    if TEShrink
        [points(iTEm1,:)]=Align3points(points([iTE,iTEm2],:),points(iTEm1,:));
        [points(iTEp1,:)]=Align3points(points([iTE,iTEp2],:),points(iTEp1,:));
    end
    
%     if points(iLEm1,2)<=points(iLE+1,2)
%          points([iLEm1,iLE+1],:)=[];
%     end
%     if points(iTEm1,2)>=points(iTEp1,2)
%          points([iTEm1,iTEp1],:)=[];
%     end
    
end

function [pointAlign]=Align3points(line,point)
    
    pointAlign=point;
    
    pointAlign(2)=(line(1,2)-line(2,2))/(line(1,1)-line(2,1))...
        *(point(1)-line(2,1))+line(2,2);
    
    if isnan(pointAlign(2)) || ~isfinite(pointAlign(2))
        pointAlign=point;
    end
    
    
end

function [indMod]=IndexMod(ind,m)
    
    indMod=mod(ind-1,m)+1;
    
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
