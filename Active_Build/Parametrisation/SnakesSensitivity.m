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
function [out]=SnakesSensitivity(entryPoint,refinedGrid,looprestart,...
        oldGrid,paramsnake,paramoptim,rootFill,varargin)
    % JUST DECLARING VARIABLES FOR LATER
    
    % plotInterval,numSteps,makeMov,boundstr
    global arrivalTolerance unstructglobal
    
    
    varExtract={'arrivalTolerance'};
    [arrivalTolerance]=ExtractVariables(varExtract,paramsnake);
    
    % ACTUALLY DOING STUFF
    refinedGriduns=ModifReshape(refinedGrid);
    
    unstructglobal=refinedGriduns;
    %profile on
    
    [out]=BuildSensitivityModes(entryPoint,refinedGriduns,refinedGrid,looprestart,...
        oldGrid,paramsnake,paramoptim,rootFill,varargin{:});
    %profile viewer
    
    
end

function [out]=...
        BuildSensitivityModes(entryPoint,refinedGriduns,refinedGrid,looprestart,...
        oldGrid,paramsnake,paramoptim,rootFill,varargin)
    % Main execution container for Snakes
    
    % Unpacking NECESSARY variables
    global snaxInitPos
    
    varExtract={'snaxInitPos'};
    [snaxInitPos]=ExtractVariables(varExtract,paramsnake);
    varExtract={'diffStepSize'};
    [diffStepSize]=ExtractVariables(varExtract,paramoptim);
    
    forceparam=paramsnake.snakes.force;
    
    
    [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        RestartSnakeProcess(looprestart);
    [cellCentredCoarse]=CellCentredSnaxelInfo(snaxel,refinedGrid,...
        cellCentredGrid,CellCentreGridInformation(oldGrid),volfracconnec);
    %CheckResultsLight(refinedGriduns,snakposition,snaxel)
    switch entryPoint
        case 'fill'
            fprintf('\n    Build Fill information from sensitivity ')
            [cellordstruct]=BuildCellConnectivity(snaxel,oldGrid,refinedGrid,volfracconnec);
            fprintf('.')
            [cellordstruct]=BuildModeMixing(cellordstruct,paramsnake);
            [cellordstruct]=ScaleModeSnakeProp(cellordstruct,cellCentredCoarse,paramsnake,diffStepSize);
            fprintf('.')
            [cellordstruct]=BuildFillMode(oldGrid,cellordstruct);
            fprintf('.')
            [newFill]=GenerateSensitivityFillPop(cellordstruct,rootFill,diffStepSize,paramoptim);
            fprintf('.')
            [out]=RemoveModeNotDesign(paramoptim,oldGrid,newFill);
            fprintf(' done!\n')
        case 'profile'
            % Eventually get gradient
            [out]=GenerateAnalyticalLoop(snaxel,refinedGriduns,refinedGrid,volfracconnec,...
                cellCentredGrid,insideContourInfo,forceparam,varargin{1},...
                oldGrid,rootFill,paramsnake,paramoptim);
    end
end


function [supportstruct]=GenerateAnalyticalLoop(snaxel,refinedGriduns,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo,forceparam,newFill,oldGrid,rootFill,paramsnake,paramoptim)
    
    varExtract={'sensAnalyticalType'};
    [sensAnalyticalType]=ExtractVariables(varExtract,paramoptim);
    
    fprintf('\n    Build Profile information from sensitivity ')
    callerStr=['GetSnaxelSensitivities(snaxel,refinedGriduns,refinedGrid,volfracconnec,',...
        'cellCentredGrid,insideContourInfo,forceparam)'];
    [~,snaxel,snakposition,sensSnax,volumefraction]=evalc(callerStr);
    %[snaxmode]=ExecuteSensitivity('analytical',snaxel,snakposition,sensSnax,volumefraction,false);
    fprintf('.')
    [snaxmode]=BuildAnalyticalMode(snaxel,snakposition,sensSnax,volumefraction,false);
    fprintf('.')
    [sensSnax]=ScaleTrimSensSnax(sensSnax,volumefraction,snaxmode,oldGrid);
    fprintf('.')
    [snaxmove]=BuildMovementStructures(snaxel,snakposition,snaxmode,sensSnax,...
        volumefraction,sensAnalyticalType);
    fprintf('.')
    deltaFill=vertcat(newFill(:).fill)-ones([numel(newFill),1])*rootFill;
    [newloop]=MoveToFill(snaxmove,deltaFill);
    fprintf('.')
    
    activeCellInd=oldGrid.cell(logical([oldGrid.cell(:).isactive])).index;
    activeSub=FindObjNum([],activeCellInd,[volumefraction(:).oldCellInd]);
    for ii=1:numel(newFill)
        supportstruct(ii).loopsens=FinishLoops(newloop(ii,:),paramsnake);
        supportstruct(ii).volumefraction=volumefraction;
        for jj=1:numel(activeSub)
            supportstruct(ii).volumefraction(activeSub(jj)).targetfill=newFill(ii).fill(jj);
        end
    end
    fprintf(' done!\n')
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

function [snaxel,snakposition,sensSnax,volumefraction]=GetSnaxelSensitivities(snaxel,refinedGriduns,refinedGrid,...
        volfracconnec,cellCentredGrid,insideContourInfo,forceparam)
    
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snakposition]=SnaxelNormal2(snaxel,snakposition);
    [volumefraction,coeffstructure,~]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo);
    forceparam.isLast=true;
%     forceparam.lengthEpsilon=0;
%     forceparam.distEpsilon=0;
    [snaxel,snakposition,~,~,sensSnax]...
        =VelocityCalculationVolumeFraction(snaxel,snakposition,volumefraction,...
        coeffstructure,forceparam);
    sensSnax(:,find(any(sensSnax~=0)))=sensSnax(:,find(any(sensSnax~=0)))./...
        repmat(sqrt(sum(sensSnax(:,find(any(sensSnax~=0))).^2)),[size(sensSnax,1),1]);
    
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


function [newGrid,newRefGrid,newRestart]=ReFillGrids(baseGrid,refinedGrid,...
        restartsnake,connectstructinfo,newFill)
    
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

function [cellordstruct]=BuildFillMode(baseGrid,cellordstruct)
    
    activeCell=logical([baseGrid.cell(:).isactive]);
    activeInd=[baseGrid.cell((activeCell)).index];
    
    
    for ii=1:length(cellordstruct)
        cellordstruct(ii).fillmode=zeros([1 numel(activeInd)]);
        sub=FindObjNum([],cellordstruct(ii).coeffMode(:,1)',activeInd);
        subAct=sub~=0;
        cellordstruct(ii).fillmode(sub(subAct))=...
            cellordstruct(ii).coeffMode(find(subAct),2)';
    end
    
    
end

function [vec]=RemoveZerosVec(vec)
    
    vec=vec(vec~=0);
    
end

function [newfillstruct]=RemoveModeNotDesign(paramoptim,baseGrid,newfillstruct)
    varExtract={'notDesInd'};
    [notDesInd]=ExtractVariables(varExtract,paramoptim);
    
    activeCell=logical([baseGrid.cell(:).isactive]);
    activeInd=[baseGrid.cell((activeCell)).index];
    kk=1;
    rmFill=[];
    for ii=1:length(newfillstruct)
        if FindObjNum([],newfillstruct(ii).cellind,activeInd(notDesInd))~=0;
            rmFill(kk)=ii;
            kk=kk+1;
        elseif FindObjNum([],newfillstruct(ii).cellind,activeInd)==0;
            rmFill(kk)=ii;
            kk=kk+1;
        end
    end
    newfillstruct(rmFill)=[];
end

function [newfillstruct]=GenerateSensitivityFillPop(cellordstruct,rootFill,diffStepSize,paramoptim)
    
    kk=1;
    
    actCell=[cellordstruct(:).index];
    actCell=RemoveIdenticalEntries(actCell);
    
    allCell=[cellordstruct(:).index];
    
    newfillstruct=repmat(struct('fill',[],'optimdat',struct('var',[],'val',[])),...
        [1,numel(diffStepSize)*numel(actCell)]);
    for ii=1:numel(diffStepSize)
        for ll=1:numel(actCell)
            jj=FindObjNum([],actCell(ll),allCell);
            fillMode=vertcat(cellordstruct(jj).fillmode);
            % fillMode=sum(fillMode,1);
            % fillMode=fillMode/max(fillMode);
            [~,pos]=max(abs(fillMode),[],1);
            fillMode=fillMode(sub2ind(size(fillMode),pos,1:size(fillMode,2)));
            newRootFill=rootFill+fillMode*diffStepSize(ii);
            [newRootFill]=OverflowHandling(paramoptim,newRootFill);
            deltaFill=newRootFill-rootFill;
            newfillstruct(kk).fill=newRootFill;
            newfillstruct(kk).optimdat.var=find(deltaFill~=0);
            newfillstruct(kk).optimdat.val=deltaFill(deltaFill~=0);
            newfillstruct(kk).cellind=actCell(ll);
            kk=kk+1;
        end
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

function [freezeVertex,borderEdgInd]=IdentifyProfileEdgeVertices(refinedGrid)
    
    borderEdges=[refinedGrid.edge(:).boundaryis0] & [refinedGrid.edge(:).boundaryis1];
    borderEdgInd=[refinedGrid.edge(borderEdges).index];
    borderVertices=sort([refinedGrid.edge(borderEdges).vertexindex]);
    freezeVertex=RemoveIdenticalEntries(borderVertices);
    
end


%% Sensitivity method development


%% Analytical snaxel Mode structure

function [snaxmode]=BuildAnalyticalMode(snaxel,snakposition,sensSnax,volumefraction,isPlot)
    
    [snaxmode]=ExtractSensitivity(snaxel,snakposition,sensSnax,volumefraction,isPlot);
    
    rootloop=OrderSurfaceSnaxel(snaxel);
    for jj=1:numel(rootloop)
        A(jj)=abs(CalculatePolyArea(rootloop(jj).snaxel.coord));
    end
    startVol=sum(A);
    
    for ii=1:numel(snaxmode)
        for jj=1:numel(snaxmode(ii).loopsnaxel)
            snaxmode(ii).A(jj)=abs(CalculatePolyArea(snaxmode(ii).loopsnaxel(jj).snaxel.coord));
        end
        [snaxmode(ii).vecmode]=BuildVectors(rootloop,snaxmode(ii).loopsnaxel);
        snaxmode(ii).deltaFrac=(sum(snaxmode(ii).A)-startVol)/snaxmode(ii).cellVol;
    end
    snaxmode([snaxmode(:).deltaFrac]==0)=[];
end

function [vectormode]=BuildVectors(looproot,loopsnaxel)
    
    for ii=1:numel(loopsnaxel)
        actRoot=FindObjNum([],loopsnaxel(ii).snaxel.index,looproot(ii).snaxel.index);
        subAct=SubDivision(loopsnaxel(ii).snaxel.coord,3,'chaikin',[0 0],'none');
        subRoot{ii}=SubDivision(looproot(ii).snaxel.coord(actRoot,:),3,'chaikin',[0 0],'none');
        vecMode{ii}=subAct-subRoot{ii};
    end
    vectormode.points=vertcat(subRoot{:});
    vectormode.vector=vertcat(vecMode{:});
end

function [snaxmode]=ExtractSensitivity(snaxel,snakposition,sensSnax,volumefraction,isPlot)
    
    if nargin<5
        isPlot=false;
    end
    global unstructglobal
    
    cellVols=[volumefraction(find(sum(abs(sensSnax))~=0)).totalvolume];
    cellInd=[volumefraction(find(sum(abs(sensSnax))~=0)).oldCellInd];
    sensSnax(:,find(sum(abs(sensSnax))==0))=[];
    dCurr=[snaxel(:).d];
    
    [snaxOrd]=SplitSnaxLoops(snaxel); % Isolates individual loops
    maxDistRatio=1/1000*min(1-abs(dCurr*2-1));
    [dChange,sensSnaxRatio]=FindModalDistanceChange(sensSnax,maxDistRatio);
    
    
    coord1=vertcat(snakposition(:).coord);
    %     [dir]=sum((vertcat(snakposition(:).vector)~=0).*[ones([length(snaxel), 1]),...
    %         ones([length(snaxel), 1])*2],2);
    [~,dir]=max(abs(vertcat(snakposition(:).vector)),[],2);
    
    for ii=1:length(sensSnax(1,:)),
        snaxCopy(ii,:)=snaxel;
        
        dAct=dCurr'+dChange{ii};
        for jj=1:length(snaxel)
            snaxCopy(ii,jj).d=dAct(jj);
        end
        [snakposition2(ii,:)]=PositionSnakes(snaxCopy(ii,:),unstructglobal);
    end
    
    coordBase{numel(snaxOrd)}=[];
    refPos{numel(snaxOrd)}=[];
    coordNew{numel(snaxOrd),numel(dChange)}=[];
    newPos{numel(snaxOrd),numel(dChange)}=[];
    newOrd{numel(snaxOrd),numel(dChange)}=[];
    for ii=1:numel(snaxOrd)
        coordBase{ii}=vertcat(snakposition(snaxOrd{ii}).coord);
        [refPos{ii}]=CreateComparisonMatrix(coordBase{ii});
        for jj=1:numel(dChange)
            coordNew{ii,jj}=vertcat(snakposition2(jj,snaxOrd{ii}).coord);
            [newPos{ii,jj}]=CreateComparisonMatrix(coordNew{ii,jj});
            [newOrd{ii,jj}]=(CompareTestpos(refPos{ii},newPos{ii,jj},1:numel(snaxOrd{ii}),dir(snaxOrd{ii})));
        end
    end
    
    for jj=1:numel(dChange)
        snaxmode(jj).cellindex=cellInd(jj);
        snaxmode(jj).snaxel=snaxCopy(jj,:);
        snaxOrdTemp=newOrd(:,jj);
        for ii=1:numel(snaxOrd)
            snaxOrdTemp{ii}=snaxOrd{ii}(snaxOrdTemp{ii});
        end
        [snaxmode(jj).snaxel]=MatchOrders(snaxmode(jj).snaxel,snaxOrdTemp);
        snaxmode(jj).loopsnaxel=OrderSurfaceSnaxel(snaxmode(jj).snaxel);
        snaxmode(jj).cellVol=cellVols(jj);
        snaxmode(jj).sensSnaxRatio=sensSnaxRatio(jj);
    end
    
    
    
    
    if isPlot
        for jj=1:numel(dChange)
            coordFull{jj}=vertcat(snakposition2(jj,:).coord);
        end
        
        for jj=1:numel(dChange)
            figure
            hold on
            for ii=1:numel(snaxOrd)
                
                plot(coordBase{ii}(:,1),coordBase{ii}(:,2),'+-')
                plot(coordNew{ii,jj}(newOrd{ii,jj},1),coordNew{ii,jj}(newOrd{ii,jj},2),'o-')
                
            end
            for ii=1:size(coord1,1)
                plot([coord1(ii,1),coordFull{jj}(ii,1)],[coord1(ii,2),coordFull{jj}(ii,2)],'k--')
            end
            
        end
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

function [snaxel]=MatchOrders(snaxel,snaxOrd)
    
    for ii=1:numel(snaxOrd)
        
        nOrd=numel(snaxOrd{ii});
        for jj=1:nOrd
            ind=mod(jj-1,nOrd)+1;
            indp1=mod(jj,nOrd)+1;
            indm1=mod(jj-2,nOrd)+1;
            
            snaxel(snaxOrd{ii}(ind)).snaxprec=snaxel(snaxOrd{ii}(indm1)).index;
            snaxel(snaxOrd{ii}(ind)).snaxnext=snaxel(snaxOrd{ii}(indp1)).index;
            snaxel(snaxOrd{ii}(ind)).connectivity=...
                [snaxel(snaxOrd{ii}(indm1)).index,snaxel(snaxOrd{ii}(indp1)).index];
        end
        
    end
    
    snaxel=snaxel([snaxOrd{:}]);
    
end

function [snaxOrd]=SplitSnaxLoops(snaxel)
    % Splits a snake into its component loops
    
    kk=1;
    jj=1;
    snaxInd=[snaxel(:).index];
    %snaxOrd=zeros(size(snaxel));
    ordList=1:numel(snaxel);
    ll=1;
    for ii=1:length(snaxel)
        kk=FindObjNum([],[snaxel(kk).snaxnext],snaxInd);
        if ordList(kk)==0
            jj=jj+1;
            kk=min(ordList(ordList~=0));
            ll=1;
        end
        snaxOrd{jj}(ll)=kk;
        ordList(kk)=0;
        ll=ll+1;
    end
    
end

function [snaxOrd]=CloseRepeatingSnaxLoops(snaxOrd)
    ii=1;
    jj=1;
    while (ii<=numel(snaxOrd))
        
        activeLoop=snaxOrd{ii};
        subActCell=sort(FindObjNum([],activeLoop(jj),activeLoop));
        if numel(subActCell)>1
            snaxOrd{ii}=activeLoop(subActCell(1):subActCell(2)-1);
            snaxOrd{numel(snaxOrd)+1}=activeLoop([subActCell(2):end,1:subActCell(1)-1]);
            jj=0;
        end
        
        ii=ii+floor(jj/length(activeLoop));
        jj=mod(jj,length(activeLoop))+1;
    end
end

function [dChange,sensSnaxRatio]=FindModalDistanceChange(sensSnax,maxDistRatio)
    
    nMode=size(sensSnax,2);
    sensSnaxRatio=zeros([1,nMode]);
    dChange{nMode}=[];
    for ii=1:nMode
        sensSnaxRatio(ii)=1/max(abs(sensSnax(:,ii)))*maxDistRatio;
        dChange{ii}=sensSnax(:,ii)*sensSnaxRatio(ii);
    end
end

function [testPos]=CreateComparisonMatrix(coord)
    
    [m,n]=size(coord);
    testPos{n}=[];
    for jj=1:n
        testPos{jj}=zeros(m);
        for ii=1:m
            testPos{jj}(:,ii)=(coord(ii,jj)>coord(:,jj))+(coord(ii,jj)<coord(:,jj))*-1;
        end
    end
    
    
end

function [newOrd]=CompareTestpos(t1,t2,ord,dir)
    
    oldOrd=ord;
    ord=FindObjNum([],ord,ord)';
    
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
    
    [newOrd,tX]=ReorderList(ord,tX);
    newOrd=newOrd(end:-1:1);
    [newOrd,tX]=ReorderList(newOrd,tX);
    newOrd=newOrd(end:-1:1);
    
    %newOrd=[newOrd(end-1),newOrd];
    rmSnak=[];
    kk=1;
    for ii=1:length(newOrd)
        
        ind=mod(ii-1,length(newOrd))+1;
        indp1=mod(ii,length(newOrd))+1;
        indm1=mod(ii-2,length(newOrd))+1;
        
        
        tTest=tDel([newOrd(indm1),newOrd(ind),newOrd(indp1)],[newOrd(ind)]);
        tTest2=t2{abs(dir(newOrd(ind))-3)}([newOrd(indm1),newOrd(ind),newOrd(indp1)]...
            ,[newOrd(ind)]);
        
        if (sum(tTest~=dir(newOrd(ind)) & tTest~=0)) ... % overtaken by neighbour
                || (sum(tTest2==0 & (tTest==dir(newOrd(ind))))) % neighbours on same line have crossed
            rmSnak(kk)=ind;
            kk=kk+1;
            
            
        end
    end
    newOrd(rmSnak)=[];
    newOrd=oldOrd(newOrd);
end

function [newOrd,tX]=ReorderList(ord,tX)
    
    newOrd=ord;
    for ii=0:2*length(ord)-1
        ind1=mod(ii,length(ord))+1;
        ind2=mod(ii+1,length(ord))+1;
        if tX(newOrd(ind2),newOrd(ind1)) % if crossing
            
            tX(newOrd(ind2),newOrd(ind1))=false; % delete the crossing flag
            tX(newOrd(ind1),newOrd(ind2))=false;
            
            interim=newOrd(ind1); % invert connection
            newOrd(ind1)=newOrd(ind2);
            newOrd(ind2)=interim;
        end
    end
    
    
end

function [loopsnaxel]=OrderSurfaceSnaxel(snaxel,snaxPositions)
    % function extracting the snaxels into their separate loops
    if nargin==1
        global unstructglobal
        snaxPositions=PositionSnakes(snaxel,unstructglobal);
    end
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

function [A]=CalculatePolyArea(points)
    
    pointsVec=points';
    pointsVec=pointsVec(:);
    %plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;
    
    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
    
end

%% Combined mode data

function [cellordstruct]=BuildCellConnectivity(snaxel,baseGrid,gridRefined,volfraconnec)
    
    edgeInd=[gridRefined.edge(:).index];
    %activeCellInd=[baseGrid.cell(logical([baseGrid.cell(:).isactive])).index];
    activeCellInd=[baseGrid.cell(:).index];
    [oldCellInd,newCellInd]=NewToOldCell(volfraconnec);
    
    for ii=1:numel(snaxel)
        snaxcellinfo(ii).index=snaxel(ii).index;
        snaxcellinfo(ii).precsnax=snaxel(ii).snaxprec;
        snaxcellinfo(ii).nextsnax=snaxel(ii).snaxnext;
        snaxcellinfo(ii).newcell=gridRefined.edge(...
            FindObjNum([],snaxel(ii).edge,edgeInd)).cellindex;
        
        snaxcellinfo(ii).oldcell=(oldCellInd(...
            FindObjNum([],snaxcellinfo(ii).newcell,newCellInd)));
    end
    
    [snaxOrd]=SplitSnaxLoops(snaxel);
    
    for ii=1:length(snaxOrd)
        cellsRaw=vertcat(snaxcellinfo(snaxOrd{ii}).oldcell);
        cellOrd{ii}=OrderRawCells(cellsRaw,activeCellInd);
    end
    delCell=false(size(cellOrd));
    for ii=1:length(cellOrd)
        delCell(ii)=isempty(cellOrd{ii});
    end
    cellOrd(delCell)=[];
    [cellOrd]=CloseRepeatingSnaxLoops(cellOrd);
    
    cellordstruct=repmat(struct('index',[],'nextcell',[],'prevcell',[],...
        'loopindex',[],'looplength',[]),[1,numel([cellOrd{:}])]);
    kk=1;
    
    for ii=1:numel(cellOrd)
        nCell=numel(cellOrd{ii});
        for jj=1:nCell
            cellordstruct(kk).index=cellOrd{ii}(jj);
            cellordstruct(kk).nextcell=cellOrd{ii}(mod(jj,nCell)+1);
            cellordstruct(kk).prevcell=cellOrd{ii}(mod(jj-2,nCell)+1);
            cellordstruct(kk).loopindex=ii;
            cellordstruct(kk).looplength=nCell;
            
            kk=kk+1;
        end
    end
    
    
end

function [cellsOrdRaw]=OrderRawCells(cellsRaw,activeCellInd)
    
    
    cellsOrdRaw=zeros([1,numel(cellsRaw)]);
    [l2Ind,l1Ind]=find((ones([2,1])*cellsRaw(1,:))==(cellsRaw(2,:)'*ones([1,2])));
    nRows=size(cellsRaw,1);
    cellsOrdRaw(1)=cellsRaw(1,l1Ind(1));
    cellsOrdRaw(2)=cellsRaw(2,l2Ind(1));
    kk=3;
    for ii=2:nRows
        l1Ind=abs(l2Ind(1)-3);
        l2Ind=find(cellsRaw(ii,l1Ind)==cellsRaw(mod(ii,nRows)+1,:));
        
        if isempty(l2Ind)
            error('Non matching cell lists')
        end
        
        cellsOrdRaw(kk)=cellsRaw(ii,l1Ind);
        cellsOrdRaw(kk+1)=cellsRaw(mod(ii,nRows)+1,l2Ind(1));
        kk=kk+2;
        
    end
    cellsOrdRaw=RemoveIdenticalConsecutivePoints(cellsOrdRaw')';
    cellsOrdRaw=cellsOrdRaw(FindObjNum([],cellsOrdRaw,activeCellInd)~=0);
    if cellsOrdRaw(end)==cellsOrdRaw(1)
        cellsOrdRaw(end)=[];
    end
end

function [oldCellInd,newCellInd]=NewToOldCell(volfraconnec)
    
    newCellInd=[volfraconnec.cell(:).newCellInd];
    oldCellInd=ones(size(newCellInd));
    kk=1;
    for ii=1:length(volfraconnec.cell)
        nNew=numel(volfraconnec.cell(ii).newCellInd);
        oldCellInd(kk:kk+nNew-1)=ones([1,nNew])*volfraconnec.cell(ii).oldCellInd;
        kk=kk+nNew;
    end
end

%% Combined mode coefficients

function [cellordstruct]=BuildModeMixing(cellordstruct,param)
    
    varExtract={'modeSmoothType','modeSmoothNum'};
    [modeSmoothType,modeSmoothNum]=ExtractVariables(varExtract,param);
    
    loopInfo=RemoveIdenticalVectors([[cellordstruct(:).loopindex]',...
        [cellordstruct(:).looplength]']);
    
    for ii=1:size(loopInfo,1)
        [coeffList{loopInfo(ii,1)}]=ReturnSmoothingCoeff(modeSmoothType,...
            loopInfo(ii,2),modeSmoothNum);
    end
    
    loopInd=[cellordstruct(:).loopindex];
    cellInd=[cellordstruct(:).index];
    
    for ii=1:length(cellordstruct)
        currLoop=loopInd==loopInd(ii);
        cellordstruct(ii).coeffMode=SingleMode(cellordstruct(currLoop),cellInd(ii),...
            coeffList{loopInd(ii)},cellInd(currLoop));
    end
    
end

function [coeffMode]=SingleMode(cellordstruct,cellNum,coeffList,cellInd)
    
    coeffMode=[cellNum,1];
    modeSub=FindObjNum([],cellNum,cellInd);
    neighbour=[cellordstruct(modeSub).prevcell,cellordstruct(modeSub).nextcell];
    for ii=2:numel(coeffList)
        coeffMode(end+1:end+2,:)=[neighbour',coeffList(ii)*ones([2,1])];
        modeSub=FindObjNum([],neighbour,cellInd);
        neighbour=[cellordstruct(modeSub(1)).prevcell,cellordstruct(modeSub(2)).nextcell];
    end
    
end

function [coeffList]=ReturnSmoothingCoeff(smoothType,nCell,nSmooth)
    % valid smoothtypes: peaksmooth polysmooth polypeaksmooth
    
    [coeffstructure]=smoothcoeffdata();
    coeffstructure=coeffstructure.(smoothType);
    nMode=numel(coeffstructure);
    for ii=1:nMode
        coeffstructure(ii).minCell=numel(coeffstructure(ii).coeff);
    end
    
    if ~strcmp(smoothType,'peaksmooth')
        minCell=[coeffstructure(:).minCell];
        maxSmooth=[coeffstructure(:).maxsmooth];
        posValid=find((minCell<=nCell) & (maxSmooth<=nSmooth));
        [~,posMaxSmooth]=max(maxSmooth(posValid));
        
        coeffstructure=coeffstructure(posValid(posMaxSmooth));
        centrePos=find(coeffstructure.coeff==1);
        
        coeffList=coeffstructure.coeff(centrePos:end);
    else
        
        centrePos=find(coeffstructure.coeff==1);
        nAct=min([nSmooth,floor((nCell-1)/2),coeffstructure.maxsmooth]);
        coeffList=coeffstructure.coeff(centrePos:centrePos+nAct);
        
    end
    
    
    
end

function [cellordstruct]=ScaleModeSnakeProp(cellordstruct,cellCentredCoarse,param,diffStepSize)
    
    try
        [modeSmoothScale]=ExtractVariables({'modeSmoothScale'},param);
    catch
        warning('Old param fed given to sensitivity calculation defaulting to scaled mode.')
        modeSmoothScale='lengthvolnormvol';
    end
    
    diffStepSize=max(abs(diffStepSize));
    [cellInd]=[cellCentredCoarse(:).index];
    [cellVol]=[cellCentredCoarse(:).volume];
    [cellLength]=[cellCentredCoarse(:).lSnax];
    [cellAct]=logical([cellCentredCoarse(:).isactive]);
    
    switch modeSmoothScale
        case 'lengthvol'
            
            for ii=1:numel(cellordstruct)
                [cellordstruct(ii).coeffMode]=ScaleModeSnake_lengthvol(...
                    cellordstruct(ii).coeffMode,cellInd,cellVol,cellLength,...
                    cellAct,diffStepSize);
            end
        case 'lengthvolnormfill' 
            % Normalizes to the original fill mode value change
            normVec= @(x) sqrt(sum(x.^2));
            for ii=1:numel(cellordstruct)
                [coeffMode]=ScaleModeSnake_lengthvol(...
                    cellordstruct(ii).coeffMode,cellInd,cellVol,cellLength,...
                    cellAct,diffStepSize);
                cellordstruct(ii).coeffMode(:,2)=coeffMode(:,2)/...
                    normVec(coeffMode(:,2))*normVec( cellordstruct(ii).coeffMode(:,2));
            end
            
        case 'lengthvolnormvol' 
            % Normalizes to the original fill mode value change
            normVec= @(x) sqrt(sum(x.^2));
            for ii=1:numel(cellordstruct)
                [coeffMode,actVol]=ScaleModeSnake_lengthvol(...
                    cellordstruct(ii).coeffMode,cellInd,cellVol,cellLength,...
                    cellAct,diffStepSize);
                cellordstruct(ii).coeffMode(:,2)=coeffMode(:,2)/...
                    normVec(coeffMode(:,2).*actVol')*normVec(...
                    cellordstruct(ii).coeffMode(:,2).*actVol');
            end  
        case 'lengthvolnormvolfill' 
            % Normalizes to the original fill mode value change
            normVec= @(x) sqrt(sum(x.^2));
            for ii=1:numel(cellordstruct)
                [coeffMode,actVol]=ScaleModeSnake_lengthvol(...
                    cellordstruct(ii).coeffMode,cellInd,cellVol,cellLength,...
                    cellAct,diffStepSize);
                cellordstruct(ii).coeffMode(:,2)=coeffMode(:,2)/...
                    normVec(coeffMode(:,2).*actVol')*normVec(...
                    cellordstruct(ii).coeffMode(:,2))*normVec(actVol');
            end
        case 'none'
            
        otherwise
            error('Invalid value for modeSmoothScale: %s ',modeSmoothScale)
    end
    
    rootMove=zeros(size(cellordstruct));
    for ii=1:numel(cellordstruct)
        rootMove(ii)=cellordstruct(ii).coeffMode(1,2);
    end
    cellOrdSub=(FindObjNum([],cellInd(cellAct),[cellordstruct(:).index]));
    cellOrdSub=cellOrdSub(cellOrdSub~=0);
    scaleRat=((cellVol(cellOrdSub).*rootMove./cellLength(cellOrdSub)));
    scaleRat=max(scaleRat(isfinite(scaleRat)));
    for ii=1:numel(cellordstruct)
        cellordstruct(ii).coeffMode(:,2)=cellordstruct(ii).coeffMode(:,2)/...
            (cellVol(cellordstruct(ii).index==cellInd)*rootMove(ii)/...
            cellLength(cellordstruct(ii).index==cellInd))*scaleRat;
    end
end

function [coeffMode,actVol]=ScaleModeSnake_lengthvol2(coeffMode,cellInd,cellVol,...
        cellLength,cellAct,dumpin)
    
    
        actSub=FindObjNum([],coeffMode(:,1),cellInd);
        actVol=cellVol(actSub);
        actLength=cellLength(actSub);
        actCellAct=cellAct(actSub);
        actCoeff=ones(size(coeffMode(:,2)));
        n=numel(actSub);
        actCoeff(2:end)=min((actLength(2:n)./actLength(max((2:n)-2,1))),1000)...
            .*(actVol(max((2:n)-2,1))./(actVol(2:n))).*actCellAct(2:n);
        actCoeff(2:2:end)=cumprod(actCoeff(2:2:end));
        actCoeff(3:2:end)=cumprod(actCoeff(3:2:end));
        coeffMode(:,2)=coeffMode(:,2).*actCoeff;
    
end


function [coeffMode,actVol,actLength]=ScaleModeSnake_lengthvol(coeffMode,cellInd,cellVol,...
        cellLength,cellAct,diffStepSize)
    
    
        actSub=FindObjNum([],coeffMode(:,1),cellInd);
        actVol=cellVol(actSub);
        actLength=cellLength(actSub);
        actCellAct=cellAct(actSub);
        actCoeff=coeffMode(:,2);
        n=numel(actSub);
        for ii=2:n
            actCoeff(ii)=(actLength(ii)./(actLength(max((ii)-2,1))+...
                sqrt(actCoeff((max((ii)-2,1)))*actVol(max((ii)-2,1))*diffStepSize)))...
            .*(actVol(max((ii)-2,1))./(actVol(ii))).*actCellAct(ii).*...
            actCoeff((max((ii)-2,1)))/coeffMode(max((ii)-2,1),2)*actCoeff(ii);
            actCoeff(isnan(actCoeff))=0;
        end
        
        
        coeffMode(:,2)=actCoeff;
    
end

%% Build modes

function [snaxmove]=BuildMovementStructures(snaxel,snakposition,snaxmode,...
        sensSnax,volumefraction,sensAnalyticalType)
    
    if numel(snaxel)~=numel(snakposition)
        error('Dimensiom mismatch between snakposition and snaxel while building movement structure')
    elseif any([snaxel(:).index]~=[snakposition(:).index])
        error('Indices mismatch between snakposition and snaxel while building movement structure')
    end
    
    [loopsnaxel,nSnax,nLoop]=CCWLoopLength(snaxel,snakposition);
    
    [snaxmove]=BuildSnaxMoveTemplate(nSnax,nSnax,nSnax);
    switch sensAnalyticalType
        case 'smooth'
            for ii=1:nLoop
                [snaxmove(ii)]=CalculateMoveData(snaxmove(ii),snakposition,sensSnax,...
                    loopsnaxel(ii).snaxel);
            end
        case 'raw'
            for ii=1:nLoop
                
                [snaxmove(ii)]=CalculateMoveDataPure(snaxmove(ii),snakposition,sensSnax,...
                    loopsnaxel(ii).snaxel);
            end
    end
    %     deltaFill=zeros([3,numel(snaxmode)]);
    %     deltaFill(1,3)=1;deltaFill(2,10)=1;deltaFill(3,25)=1;
    %     deltaFill=deltaFill*1e-2;
    %     MoveToFill(snaxmove,deltaFill)
end

function [loopsnaxel,nSnax,nLoop]=CCWLoopLength(snaxel,snakposition)
    
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel,snakposition);
    nLoop=numel(loopsnaxel);
    nSnax=zeros([1,nLoop]);
    for ii=1:nLoop
        if ~CCWLoop(loopsnaxel(ii).snaxel.coord);
            loopsnaxel(ii).snaxel.index=flip(loopsnaxel(ii).snaxel.index);
            loopsnaxel(ii).snaxel.coord=flip(loopsnaxel(ii).snaxel.coord);
        end
        % Add the length coordinates
        loopsnaxel(ii).snaxel.edgelength=sqrt(sum((loopsnaxel(ii).snaxel.coord([2:end,1],:)-...
            loopsnaxel(ii).snaxel.coord).^2,2));
        loopsnaxel(ii).snaxel.lpos=cumsum([0;loopsnaxel(ii).snaxel.edgelength]);
        nSnax(ii)=numel(loopsnaxel(ii).snaxel.index);
    end
end

function [snaxmove]=BuildSnaxMoveTemplate(nSnax,nEdge,nVertex)
    
    nLoop=numel(nSnax);
    
    snaxmove=struct('snax',struct([]),'vertex',struct([]),'edge',struct([]),...
        'support',struct('nSnax',[],'maxL',[]));
    snaxmove=repmat(snaxmove,[1,nLoop]);
    
    snaxTemp=struct('index',[],'rN',[],'rT',[],'d',[],'sens',[],'posL',[]);
    edgeTemp=struct('index',0,'posL',[],'coord',[],'normal',[]);
    vertTemp=edgeTemp;
    
    for ii=1:nLoop
        
        snaxmove(ii).snax=repmat(snaxTemp,[1,nSnax(ii)]);
        snaxmove(ii).vertex=repmat(edgeTemp,[1,nVertex(ii)]);
        snaxmove(ii).edge=repmat(vertTemp,[1,nEdge(ii)]);
        snaxmove(ii).support.nSnax=nSnax(ii);
    end
    
end

function [snaxmove]=CalculateMoveData(snaxmove,snakposition,sensSnax,loopsnaxel)
    
    fullSnax=[snakposition(:).index];
    
    snaxSub=FindObjNum([],loopsnaxel.index,fullSnax);
    snaxSub=reshape(snaxSub,[1,numel(snaxSub)]);
    
    snaxmove.support.maxL=loopsnaxel.lpos(end);
    
    
    kk=1;
    for ii=snaxSub
        kks=[kk:min([kk+1,snaxmove.support.nSnax]),...
            mod(kk,snaxmove.support.nSnax)+1:1];
        snaxmove.snax(kk).index=snakposition(ii).index;
        snaxmove.snax(kk).d=sqrt(sum(snakposition(ii).vectornotnorm.^2));
        snaxmove.snax(kk).sens=sensSnax(ii,:);
        snaxmove.snax(kk).posL=loopsnaxel.lpos(kk);
        
        snaxmove.vertex(kk).index=snakposition(ii).index;
        snaxmove.vertex(kk).normal=sum(snakposition(ii).normvector);
        snaxmove.vertex(kk).normal=snaxmove.vertex(kk).normal/...
            sqrt(sum(snaxmove.vertex(kk).normal.^2));
        snaxmove.vertex(kk).coord=snakposition(ii).coord;
        snaxmove.vertex(kk).posL=loopsnaxel.lpos(kk);
        
        snaxmove.edge(kk).index=loopsnaxel.index(kks);
        snaxmove.edge(kk).coord=mean(loopsnaxel.coord(kks,:));
        snaxmove.edge(kk).posL=mean(loopsnaxel.lpos(kk:kk+1));
        snaxmove.edge(kk).normal=snakposition(ii).vectornext;
        
        [vecAngles]=ExtractAnglepm180(snaxmove.vertex(kk).normal,snakposition(ii).vector);
        snaxmove.snax(kk).rT=sin(vecAngles);
        snaxmove.snax(kk).rN=cos(vecAngles);
        kk=kk+1;
    end
    
    
    
    
end

function [snaxmove]=CalculateMoveDataPure(snaxmove,snakposition,sensSnax,loopsnaxel)
    
    fullSnax=[snakposition(:).index];
    
    snaxSub=FindObjNum([],loopsnaxel.index,fullSnax);
    snaxSub=reshape(snaxSub,[1,numel(snaxSub)]);
    
    snaxmove.support.maxL=loopsnaxel.lpos(end);
    
    
    kk=1;
    for ii=snaxSub
        kks=[kk:min([kk+1,snaxmove.support.nSnax]),...
            mod(kk,snaxmove.support.nSnax)+1:1];
        snaxmove.snax(kk).index=snakposition(ii).index;
        snaxmove.snax(kk).d=sqrt(sum(snakposition(ii).vectornotnorm.^2));
        snaxmove.snax(kk).sens=sensSnax(ii,:);
        snaxmove.snax(kk).posL=loopsnaxel.lpos(kk);
        
        snaxmove.vertex(kk).index=snakposition(ii).index;
        snaxmove.vertex(kk).normal=snakposition(ii).vector;
        snaxmove.vertex(kk).normal=snaxmove.vertex(kk).normal/...
            sqrt(sum(snaxmove.vertex(kk).normal.^2));
        snaxmove.vertex(kk).coord=snakposition(ii).coord;
        snaxmove.vertex(kk).posL=loopsnaxel.lpos(kk);
        
        snaxmove.edge(kk).index=loopsnaxel.index(kks);
        snaxmove.edge(kk).coord=mean(loopsnaxel.coord(kks,:));
        snaxmove.edge(kk).posL=mean(loopsnaxel.lpos(kk:kk+1));
        snaxmove.edge(kk).normal=snakposition(ii).vectornext;
        
        [vecAngles]=ExtractAnglepm180(snaxmove.vertex(kk).normal,snakposition(ii).vector);
        snaxmove.snax(kk).rT=0;
        snaxmove.snax(kk).rN=1;
        kk=kk+1;
    end
    
    
    
    
end

function [sensSnax]=ScaleTrimSensSnax(sensSnax,volumefraction,snaxmode,baseGrid)
    % Scale the sensitivity values in sensSnax to be in units of d per unit
    % of volume fraction
    
    
    cellInd=[volumefraction(:).oldCellInd];
    cellSens=[snaxmode(:).cellindex];
    
    sensRatio=[snaxmode(:).sensSnaxRatio]./[snaxmode(:).deltaFrac];
    colSub=FindObjNum([],cellSens,cellInd);
    
    sensSnax(:,colSub)=sensSnax(:,colSub).*(ones([size(sensSnax,1),1])*sensRatio);
    
    activeCell=logical([baseGrid.cell(:).isactive]);
    activeSub=find(activeCell);
    
    sensSnax=sensSnax(:,activeSub);
end

% Movement function

function [newloop]=MoveToFill(snaxmove,deltaFill)
    % deltaFill is a matrix where each row is a different fill
    
    
    % Find new Lpos
    
    nFill=size(deltaFill,1);
    newloop=repmat(struct('snaxel',struct('index',[],'coord',[],'vector',[],...
        'coordnoscale',[])),[nFill,numel(snaxmove)]);
    
    % deltaFill must be trimmed to only have active snax boxes
    for ii=1:numel(snaxmove)
        % each column represent a different fill
        % each row represents a different snaxel
        
        % preparation
        nSnax=numel(snaxmove(ii).snax);
        nEdge=numel(snaxmove(ii).edge);
        nVert=numel(snaxmove(ii).vertex);
        snaxList=[snaxmove(ii).snax(:).index];
        vertPosL=repmat(reshape([snaxmove(ii).vertex(:).posL],...
            [1 1 nVert]),[nSnax,nFill,1]);
        edgePosL=repmat(reshape([snaxmove(ii).edge(:).posL],...
            [1 1 nEdge]),[nSnax,nFill,1]);
        
        % vertex recalc
        mvtCoeffFill=(vertcat(snaxmove(ii).snax(:).sens)*deltaFill')...
            .*(vertcat(snaxmove(ii).snax(:).d)*(ones([1,size(deltaFill,1)])));
        
        coeffNormpos=mvtCoeffFill.*(vertcat(snaxmove(ii).snax(:).rN)*(ones([1,size(deltaFill,1)])));
        
        deltaLpos=mvtCoeffFill.*(vertcat(snaxmove(ii).snax(:).rT)*(ones([1,size(deltaFill,1)])));
        newLpos=deltaLpos+(vertcat(snaxmove(ii).snax(:).posL)*(ones([1,size(deltaFill,1)])));
        newLpos(newLpos<0)=snaxmove(ii).support.maxL+newLpos(newLpos<0);
        newLpos(newLpos>snaxmove(ii).support.maxL)=...
            newLpos(newLpos>snaxmove(ii).support.maxL)-snaxmove(ii).support.maxL;
        [~,newOrder]=sort(newLpos);
%         plot(newOrder)
        [distVert,closeVert]=min(abs(vertPosL-repmat(newLpos,[1 1 nVert])),[],3);
        [distEdge,closeEdge]=min(abs(edgePosL-repmat(newLpos,[1 1 nEdge])),[],3);
        
        for jj=1:nFill
            
            vertCoord=reshape(vertcat(snaxmove(ii).vertex(closeVert(:,jj)).coord),[nVert,2,1]);
            edgeCoord=reshape(vertcat(snaxmove(ii).edge(closeEdge(:,jj)).coord),[nVert,2,1]);
            newCoordTemp=(vertCoord.*(distEdge(:,jj)*ones([1 2]))...
                +edgeCoord.*(distVert(:,jj)*ones([1 2])))...
                ./((distVert(:,jj)+distEdge(:,jj))*ones([1 2]));
            newLCoord=newCoordTemp;
            
            vertNormal=reshape(vertcat(snaxmove(ii).vertex(closeVert(:,jj)).normal),[nVert,2,1]);
            edgeNormal=reshape(vertcat(snaxmove(ii).edge(closeEdge(:,jj)).normal),[nVert,2,1]);
            newNormalTemp=(vertNormal.*(distEdge(:,jj)*ones([1 2]))...
                +edgeNormal.*(distVert(:,jj)*ones([1 2])))...
                ./((distVert(:,jj)+distEdge(:,jj))*ones([1 2]));
            newNormal=newNormalTemp./(sqrt(sum(newNormalTemp.^2,2))*ones([1 2]));
            
            newCoord=(coeffNormpos(:,jj)*ones([1 2])).*...
                newNormal+newLCoord;
            
            newloop(jj,ii).snaxel.index=snaxList(newOrder(:,jj));
            newloop(jj,ii).snaxel.coord=newCoord(newOrder(:,jj),:);
            newloop(jj,ii).snaxel.coordnoscale=newCoord(newOrder(:,jj),:);
            newloop(jj,ii).snaxel.vector=newNormal(newOrder(:,jj),:);
        end
        
    end
    
    
    if false
        
        plotPoints=@(points,str) plot(points([1:end,1],1),points([1:end,1],2),str);
        quiverPoints= @(points,vector,str) quiver(points([1:end,1],1),...
            points([1:end,1],2),vector([1:end,1],1),vector([1:end,1],2));
        figure,hold on
        axis equal
        
        
        colorStr='rky';
        
        for jj=1:nFill
            for ii=1:numel(snaxmove),
                plotPoints(newloop(jj,ii).snaxel.coord,[colorStr(jj),'-s'])
            end
        end
        for ii=1:numel(snaxmove),
            plotPoints(vertcat(snaxmove(ii).vertex(:).coord),'b-+')
            quiverPoints(vertcat(snaxmove(ii).vertex(:).coord),vertcat(snaxmove(ii).vertex(:).normal))
            plotPoints(vertcat(snaxmove(ii).edge(:).coord),'go')
            quiverPoints(vertcat(snaxmove(ii).edge(:).coord),vertcat(snaxmove(ii).edge(:).normal))
        end
    end
end