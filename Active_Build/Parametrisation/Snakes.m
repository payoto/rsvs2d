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
function [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=Snakes(refinedGrid,looprestart,...
        oldGrid,connectionInfo,param)
    % JUST DECLARING VARIABLES FOR LATER
    
    % plotInterval,numSteps,makeMov,boundstr
    global arrivalTolerance unstructglobal maxStep maxDt
    
    
    varExtract={'arrivalTolerance','maxStep','maxDt','snakesConsole','case'};
    [arrivalTolerance,maxStep,maxDt,snakesConsole,caseStr]=ExtractVariables(varExtract,param);
    
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
        figure('Name',['VolError',caseStr]),semilogy(1:length(snakSave),[snakSave(:).currentConvVolume])
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
    
    varExtract={'mergeTopo','boundstr','snakesConsole','dtRatio','snaxInitPos','checkSensitivities','refineGrid'};
    [mergeTopo,boundstr,snakesConsole,dtRatio,snaxInitPos,checkSensitivities,refineGrid]=ExtractVariables(varExtract,param);
    
    forceparam=param.snakes.force;
    dtMin=maxDt/dtRatio;
    
    
    % Starting process
    disp(['  Start initialisation'])
    tStepStart=now;
    
    
    edgeOrient=[refinedGrid.edge(:).orientation];
    edgeIndList=[refinedGrid.edge(:).index];
    if numel(refineGrid)==1
        refGridRatio=1;
    else
        refGridRatio=refineGrid(1)/refineGrid(2);
    end
    edgeDat.orient=edgeOrient;
    edgeDat.indList=edgeIndList;
    edgeDat.ratio=refGridRatio;
    
    
    [cellCentredGrid,volfracconnec,borderVertices,snaxel,...
        insideContourInfo]=InitialisationRestart(refinedGriduns,...
        refinedGrid,looprestart,oldGrid,connectionInfo,mergeTopo,boundstr,param,edgeDat);
    
    
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snakposition]=SnaxelNormal2(snaxel,snakposition);
    %CheckResultsLight(refinedGriduns,snakposition,snaxel)
    
    
    [freezeVertex,borderEdgInd]=IdentifyProfileEdgeVertices(refinedGrid);
    
    tStepEnd=now;
    disp(['  Initialisation time:',datestr(tStepEnd-tStepStart,'HH:MM:SS:FFF')]);
    tStart=now;
    
    disp(['  Start Time Stepping'])
    
    if snakesConsole
        [ii,snaxel,snakposition,insideContourInfo,forceparam,snakSave,currentConvVolume,edgeDat]=...
            IterSnakes(param,snaxel,refinedGriduns,refinedGrid,volfracconnec,cellCentredGrid,...
            insideContourInfo,forceparam,oldGrid,maxStep,maxDt,dtMin,borderVertices);
    else
        [~,ii,snaxel,snakposition,insideContourInfo,forceparam,snakSave,currentConvVolume,edgeDat]...
            =evalc(['IterSnakes(param,snaxel,refinedGriduns,refinedGrid,volfracconnec,',...
            'cellCentredGrid,insideContourInfo,forceparam,oldGrid,maxStep,maxDt,dtMin,borderVertices);']);
    end
    
    tEnd=now;
    disp(['  ',int2str(ii),' Steps Performed'])
    disp(['  Iteration Time:',datestr(tEnd-tStart,'HH:MM:SS:FFF')]);
    disp(['  Volume converged to ',num2str(currentConvVolume,'%.5e')])
    
    [snaxel,snakposition,loopsnaxel]=FinishSnakes(snaxel,...
        borderVertices,refinedGriduns,param,edgeDat);
    
    if checkSensitivities
        GetSnaxelSensitivities(snaxel,refinedGriduns,refinedGrid,volfracconnec,...
            cellCentredGrid,insideContourInfo,forceparam);
    end
    if snakesConsole
        %CheckResultsLight(refinedGriduns,snakposition,snaxel)
    end
    restartsnake.snaxel=snaxel;
    restartsnake.insideContourInfo=insideContourInfo;
    restartsnake.cellCentredGrid=cellCentredGrid;
    restartsnake.volfracconnec=volfracconnec;
    restartsnake.borderVertices=borderVertices;
    
    
    
end

function [ii,snaxel,snakposition,insideContourInfo,forceparam,snakSave,currentConvVolume,edgeDat]=...
        IterSnakes(param,snaxel,refinedGriduns,refinedGrid,volfracconnec,cellCentredGrid,...
        insideContourInfo,forceparam,oldGrid,maxStep,maxDt,dtMin,borderVertices)
    
    global snaxInitPos
    varExtract={'snakesSteps','mergeTopo','makeMov','convLevel','debugPlot','plotInterval',...
        'subStep','snakesMinSteps','stepType','vSwitch','convCheckRate',...
        'convCheckRange','convDistance','dtRatio','refineGrid'};
    [snakesSteps,mergeTopo,makeMov,convLevel,debugPlot,...
        plotInterval,subStep,snakesMinSteps,stepType,vSwitch,convCheckRate,...
        convCheckRange,convDistance,dtRatio,refineGrid]...
        =ExtractVariables(varExtract,param);
    
    trigCount=0;
    movFrame=struct([]);
    isConvergingPast=true;
    isChangeddtMax=false;
    lastConvCheck=0;
    dtMinStart=dtMin;
    [nonBreedVertPersist]=SetNonBreedVerticesPersistent(param,refinedGrid);
    borderVertices.nonbreedpers=nonBreedVertPersist;
    edgeOrient=[refinedGrid.edge(:).orientation];
    edgeIndList=[refinedGrid.edge(:).index];
    if numel(refineGrid)==1
        refGridRatio=1;
    else
        refGridRatio=refineGrid(1)/refineGrid(2);
    end
    edgeDat.orient=edgeOrient;
    edgeDat.indList=edgeIndList;
    edgeDat.ratio=refGridRatio;
    cntSave=1;
    for ii=1:snakesSteps
        %snaxInitPos=snaxInitPos*min(exp(-1/20*(ii-snakesSteps/2)),1);
        
        fprintf('     Step %i  -',ii);
        tStepStart=now;
        %snaxel=SnaxelDistanceUpdate(snaxel,0.1,ones([1,length(snaxel)]),ones([1,length(snaxel)]));
        %arrivalTolerance=arrivalTolerance1*exp(-decayCoeff*ii);
        
        % snaxel properties calculation
        [snakposition]=PositionSnakes(snaxel,refinedGriduns);
        [snakposition]=SnaxelNormal2(snaxel,snakposition);
        
        [volumefraction,coeffstructure,cellCentredGridSnax,convergenceCondition,...
            currentConvVelocity,currentConvVolume,forceparam,lastAlgo,trigCount,...
            snaxel,snakposition,snaxelmodvel,velcalcinfo]=...
            VelocityAndVolumeProcess(param,snaxel,snakposition,refinedGrid,volfracconnec,...
            cellCentredGrid,insideContourInfo,convLevel,forceparam,ii,trigCount);
        
        
        
        if convergenceCondition && ii>snakesMinSteps && lastAlgo
            
            fprintf(' -  Snakes Converged!\n')
            break
        end
        fprintf([' DVol = ', num2str(currentConvVolume)]);
        % visualise results
        if ((round(ii/plotInterval)==ii/plotInterval) && plotInterval) || sum(ii==debugPlot)
            [movFrame]=CheckResults(ii,refinedGriduns,oldGrid,snakposition,snaxelmodvel,makeMov,volumefraction);
            %[movFrame]=CheckResults(ii,refinedGriduns,oldGrid,snakposition,snaxel,makeMov,volumefraction);
            
        end
        
        [dt,dtSnax,maxDist]=TimeStepCalculation(snaxel,maxStep,maxDt,dtMin,...
            stepType,vSwitch,edgeDat,mergeTopo);
        % Save and exit conditions
        [snakSave(cntSave)]=WriteSnakSave(param,snaxel,dt,snakposition,...
            volumefraction,cellCentredGridSnax,currentConvVelocity,...
            currentConvVolume,movFrame,velcalcinfo,insideContourInfo);
        cntSave=cntSave+1;
        
        if mod(ii,convCheckRate)==0 && ii>convCheckRange ...
                && (~isChangeddtMax || (mod(ii-lastConvCheck,convCheckRate)==0))
            [maxDt,isConvergingPast,isChangeddtMax]=TestConvergenceRate(...
                [snakSave((end-convCheckRange):end).currentConvVolume],convLevel,maxDt,...
                convDistance,isConvergingPast);
            maxDt=max([maxDt,dtMinStart]);
            dtMin=maxDt/dtRatio;
            lastConvCheck=ii;
        end
        
        for subStep=1:subStep
            snaxel=SnaxelDistanceUpdate(snaxel,dt,dtSnax,maxDist,stepType);
            [snakposition]=PositionSnakes(snaxel,refinedGriduns);
            %             [snakposition]=SnaxelNormal2(snaxel,snakposition);
            %             [volumefraction,coeffstructure,cellCentredGridSnax,convergenceCondition,...
            %                 currentConvVelocity,currentConvVolume,forceparam,lastAlgo,trigCount,...
            %                 snaxel,snakposition,snaxelmodvel,velcalcinfo]=...
            %                 VelocityAndVolumeProcess(snaxel,snakposition,refinedGrid,volfracconnec,...
            %                 cellCentredGrid,insideContourInfo,convLevel,forceparam,ii,trigCount);
        end
        % Topology Trimming, Merging and Freezing
        [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
        [snakposition]=PositionSnakes(snaxel,refinedGriduns);
        [snaxel]=FreezingFunction(snaxel,borderVertices,edgeDat,mergeTopo);
        [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo);
        
%         [snakSave(cntSave)]=WriteSnakSave(param,snaxel,dt,snakposition,...
%             volumefraction,cellCentredGridSnax,currentConvVelocity,...
%             currentConvVolume,movFrame,velcalcinfo,insideContourInfo);
%         cntSave=cntSave+1;
        % Snaxel Repopulation In both directions
        
        [nonBreedVert]=SetNonBreedVertices(borderVertices,ii,param);
        nonBreedVert=[nonBreedVert,nonBreedVertPersist];
        borderVertices.nonBreedVert=nonBreedVert;
        [snaxel,insideContourInfo,nonBreedVert]=SnaxelBreeding(snaxel,...
            insideContourInfo,refinedGriduns,nonBreedVert,edgeDat,mergeTopo);
        if numel(snaxel)==0
            warning('Contour has collapsed')
            break
        end
        
        [snaxelrev,insideContourInfoRev]=ReverseSnaxelInformation(snaxel,...
            insideContourInfo,refinedGriduns);
        
        [snaxelrev,insideContourInfoRev,~]=SnaxelBreeding(snaxelrev,...
            insideContourInfoRev,refinedGriduns,nonBreedVert,edgeDat,mergeTopo);
        if numel(snaxelrev)==0
            warning('Contour has collapsed')
            break
        end
        
        [snaxel,insideContourInfo]=ReverseSnaxelInformation(snaxelrev,...
            insideContourInfoRev,refinedGriduns);
        
        
        % Topology Trimming, Merging and Freezing
        [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
        [snakposition]=PositionSnakes(snaxel,refinedGriduns);
        [snaxel]=FreezingFunction(snaxel,borderVertices,edgeDat,mergeTopo);
        [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo);
        
        
        nSnax=length(snaxel);
        frozen=sum(logical([snaxel(:).isfreeze]));
        tStepEnd=now;
        
        fprintf(['  Time Taken: ',datestr(tStepEnd-tStepStart,'HH:MM:SS:FFF'),'\n']);
        
        if nSnax==frozen
            break
        end
        
    end
    
end

function [nonBreedVert]=SetNonBreedVertices(borderVertices,iterNum,param)
    varExtract={'vertLooseStep','pinnedVertex'};
    [vertLooseStep,pinnedVertex]=ExtractVariables(varExtract,param);
    nonBreedVert=[];
    
    if iterNum<=vertLooseStep
        nonBreedVert=[nonBreedVert,borderVertices.weak,borderVertices.loop];
    end
    nonBreedVert=unique(nonBreedVert);
end

function [nonBreedVertPersist]=SetNonBreedVerticesPersistent(param,refinedGrid)
    normVec=@(vec,ref) sqrt(sum((vec-(ones(size(vec,1),1)*ref)).^2,2));
    
    varExtract={'pinnedVertex'};
    [pinnedVertex]=ExtractVariables(varExtract,param);
    nonBreedVertPersist=[];
    
    switch pinnedVertex
        case 'none'
            nonBreedVertPersist=[];
            
        case 'LETE'
            coord=vertcat(refinedGrid.vertex(:).coord);
            [~,iLE]=min(normVec(coord,[0 0]));
            [~,iTE]=min(normVec(coord,[1 0]));
            nonBreedVertPersist=[refinedGrid.vertex(iLE).index,refinedGrid.vertex(iTE).index];
    end
    
    nonBreedVertPersist=unique(nonBreedVertPersist);
end

function [cellCentredGrid,volfracconnec,borderVertices,snaxel,...
        insideContourInfo]=InitialisationRestart(refinedGriduns,...
        refinedGrid,looprestart,oldGrid,connectionInfo,mergeTopo,boundstr,param,edgeDat)
    
    varExtract={'restart'};
    [restart]=ExtractVariables(varExtract,param);
    if ~restart
        [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
            StartSnakeProcess(refinedGriduns,refinedGrid,looprestart,...
            oldGrid,connectionInfo,mergeTopo,boundstr,edgeDat);
    else
        [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
            RestartSnakeProcess(looprestart);
    end
    
    
end

function [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        StartSnakeProcess(refinedGriduns,refinedGrid,loop,...
        oldGrid,connectionInfo,mergeTopo,boundstr,edgeDat)
    disp('    Generate Cell Centred Grid')
    [cellCentredGrid]=CellCentreGridInformation(refinedGrid);
    disp('    Establish Cell Volume Information')
    [volfracconnec]=VolumeFractionConnectivity(oldGrid,...
        connectionInfo,cellCentredGrid,refinedGrid);
    
    
    insideContourInfo=refinedGriduns.edge.(boundstr{2});
    disp('    Find Border Vertices')
    [borderVertices.strong]=FindBorderVertex(refinedGriduns);
    disp('    Initialise Snaxel Grid')
    [borderVertices.weak]=FindWeakBorderVertex(oldGrid,cellCentredGrid,volfracconnec);
    % Inside contour info will change depending on the type of contour under
    % consideration (0 contour or 1 contour)
    borderVertices.loop=[];
    for ii=1:numel(loop)
        borderVertices.loop=[borderVertices.loop,reshape(loop(ii).vertex.index,[1,numel(loop(ii).vertex.index)])];
    end
    [snaxel,insideContourInfo]=SnaxelInitialisation(refinedGriduns,loop,insideContourInfo,boundstr);
    [snaxel]=InitisaliseSnaxelOrder(snaxel);
    [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snaxel,insideContourInfo]=TopologyMergingProcess(snaxel,snakposition,insideContourInfo);
    [snaxel]=FreezingFunction(snaxel,borderVertices,edgeDat,mergeTopo);
end

function [cellCentredGrid,volfracconnec,borderVertices,snaxel,insideContourInfo]=...
        RestartSnakeProcess(restartsnake)
    
    [cellCentredGrid]=restartsnake.cellCentredGrid;
    [volfracconnec]=restartsnake.volfracconnec;
    [borderVertices]=restartsnake.borderVertices;
    
    snaxel=restartsnake.snaxel;
    insideContourInfo=restartsnake.insideContourInfo;
    [snaxel]=AddOrderEdgeField(snaxel);
end

function [snaxel,snakposition,loopsnaxel]=FinishSnakes(snaxel,...
        borderVertices,refinedGriduns,param,edgeDat)
    
    
    %disp('Finished Iterations , starting Post Process')
    [snaxel]=FreezingFunction(snaxel,borderVertices,edgeDat);
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    
    %disp('Creating Snaxel Loops')
    [loopsnaxel]=ExtractSnaxelLoops(snaxel,param);
    
end


function []=GetSnaxelSensitivities(snaxel,refinedGriduns,refinedGrid,...
        volfracconnec,cellCentredGrid,insideContourInfo,forceparam)
    
    [snakposition]=PositionSnakes(snaxel,refinedGriduns);
    [snakposition]=SnaxelNormal2(snaxel,snakposition);
    [volumefraction,coeffstructure,cellCentredGridSnax]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo);
    forceparam.isLast=true;
    forceparam.lengthEpsilon=0;
    [snaxel,snakposition,snaxelmodvel,velcalcinfo,forceparam,sensSnax]...
        =VelocityCalculationVolumeFraction(snaxel,snakposition,volumefraction,...
        coeffstructure,forceparam);
    testSensitivity(snaxel,snakposition,sensSnax)
end

function [snaxelRev,insideContourInfoRev]=ReverseSnaxelInformation(snaxel,...
        insideContourInfo,unstructured)
    % Reverses all the snaxel information
    [snaxelRev]=ReverseSnakes(snaxel);
    [insideContourInfoRev]=ReverseInsideContourInfo(snaxel,...
        insideContourInfo,unstructured);
    
end

%% Iteration sub parts

function [snaxel,insideContourInfo,nonBreedVert]=SnaxelBreeding(snaxel,...
        insideContourInfo,refinedGriduns,nonBreedVert,edgeDat,mergeTopo)
    
    
    [snaxel,insideContourInfo,nonBreedVert]=SnaxelRepopulate(refinedGriduns,snaxel,...
        insideContourInfo,nonBreedVert,edgeDat,mergeTopo);
    [snaxel,insideContourInfo]=SnaxelCleaningProcess(snaxel,insideContourInfo);
end

function [convergenceCondition,currentConvVelocity,currentConvVolume]=...
        ConvergenceTest(snaxel,volumefraction,convLevel)
    % Calculates wether the snaxel population has converged
    
    % Get profile lengths
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel);
    
    
    
    % Velocity condition
    vSnax=[snaxel(~logical([snaxel(:).isfreeze])).v];
    
    currentConvVelocity=(sqrt(sum(vSnax.^2))/length(vSnax));
    conditionVelocity=currentConvVelocity<convLevel;
    % Volume fraction condition
    isCellSnax=[volumefraction(:).isSnax];
    diffVolFrac=[volumefraction(isCellSnax).targetfill]-[volumefraction(isCellSnax).volumefraction];
    currentConvVolume=(sqrt(sum(diffVolFrac.^2))/length(diffVolFrac));
    conditionVolume=currentConvVolume<convLevel;
    
    convergenceCondition= conditionVolume && conditionVelocity;
    %convergenceCondition= conditionVolume;
end

function [forceparam,lastAlgo,trigCount]=CheckCurrentVelocityAlgorithm(forceparam,currStep,convLevel,trigCount)
    
    maxTypeSub=length(forceparam.vel.Type);
    switch forceparam.vel.ChangeTrigger
        case 'none'
            typSub=1;
        case 'step'
            
            typSub=find(forceparam.vel.ChangeStep<currStep,1,'last');
        case 'conv'
            if trigCount<length(forceparam.vel.ChangeConv)
                if convLevel<forceparam.vel.ChangeConv(trigCount+1)
                    trigCount=trigCount+1;
                end
            end
            typSub=trigCount;
            
        case 'both'
            typSubStep=find(forceparam.vel.ChangeStep<currStep,1,'last');
            if trigCount<length(forceparam.vel.ChangeConv)
                if convLevel<forceparam.vel.ChangeConv(trigCount+1)
                    trigCount=trigCount+1;
                end
            end
            typSubConv=trigCount;
            typSub=max([typSubStep,typSubConv]);
    end
    
    forceparam.velType=forceparam.vel.Type{typSub};
    lastAlgo=typSub==maxTypeSub;
end

function [volumefraction,coeffstructure,cellCentredGridSnax,convergenceCondition,...
        currentConvVelocity,currentConvVolume,forceparam,lastAlgo,trigCount,...
        snaxel,snakposition,snaxelmodvel,velcalcinfo]=...
        VelocityAndVolumeProcess(param,snaxel,snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo,convLevel,forceparam,ii,trigCount)
    
    [volfracconnec,cellCentredGrid]=ModifyFillInformation...
        (cellCentredGrid,volfracconnec,param,ii-1);
    
    [volumefraction,coeffstructure,cellCentredGridSnax]=VolumeFraction(snaxel,...
        snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo);
    
    [convergenceCondition,currentConvVelocity,currentConvVolume]=...
        ConvergenceTest(snaxel,volumefraction,convLevel);
    
    [forceparam,lastAlgo,trigCount]=CheckCurrentVelocityAlgorithm(forceparam,...
        ii,currentConvVolume,trigCount);
    
    [volumefraction]=ModifyFillInformationError...
        (volumefraction,param,ii-1);
    
    [snaxel,snakposition,snaxelmodvel,velcalcinfo,forceparam]=...
        VelocityCalculationVolumeFraction(snaxel,snakposition,volumefraction,...
        coeffstructure,forceparam);
    
end

function [volfraconnec,cellCentredGrid]=ModifyFillInformation...
        (cellCentredGrid,volfraconnec,param,nStep)
    
    varExtract={'fillLooseStep','fillLooseCut'};
    [fillLooseStep,fillLooseCut]=ExtractVariables(varExtract,param);
    
    newInd=[cellCentredGrid(:).index];
    fill=[volfraconnec.cell(:).targetfill];
    if nStep<fillLooseStep
        warning(['Attempting to use deprecated option, Use fillErrStep and fillErrCut instead\n ',...
            'code needs to be uncommented if this option needs to be used'])
%         modifFillPos=find(fill<fillLooseCut & fill~=0);
%         for ii=modifFillPos
%             fillDelta=(fillLooseCut-volfraconnec.cell(ii).targetfill)-...
%                 nStep*(fillLooseCut-volfraconnec.cell(ii).targetfill)/fillLooseStep;
%             volfraconnec.cell(ii).targetfill=fillDelta+volfraconnec.cell(ii).targetfill;
%             
%             currSub=FindObjNum([],volfraconnec.cell(ii).newCellInd,newInd);
%             [cellCentredGrid(currSub).fill]=deal(fillDelta+volfraconnec.cell(ii).targetfill);
%         end
%         
%         modifFillPos=find(fill<=0);
%         for ii=modifFillPos
%             fillDelta=-0.1;
%             volfraconnec.cell(ii).targetfill=fillDelta;
%             
%             currSub=FindObjNum([],volfraconnec.cell(ii).newCellInd,newInd);
%             [cellCentredGrid(currSub).fill]=deal(fillDelta);
%         end
%         
%         modifFillPos=find(fill>=1);
%         for ii=modifFillPos
%             fillDelta=1.1;
%             volfraconnec.cell(ii).targetfill=fillDelta;
%             
%             currSub=FindObjNum([],volfraconnec.cell(ii).newCellInd,newInd);
%             [cellCentredGrid(currSub).fill]=deal(fillDelta);
%         end
    end
end


function [volumefraction]=ModifyFillInformationError(volumefraction,param,nStep)
   % visualisation of the conditions 
 conds=@(a,c) {a.*c+a,c.*(1-a)+a,-a.*c+a,a-c.*(1-a)};
% [a,c]=meshgrid(linspace(0,1,100),1:100);
% zVals=conds(a,c);
% figure,hold on
% surf(a,c,min(min(zVals{1},zVals{2}),1),'facecolor',[0 1 0])
% surf(a,c,max(max(zVals{3},zVals{4}),0),'facecolor',[0 0 1])
    
    
    varExtract={'fillErrStep','fillErrCut'};
    [fillErrStep,fillErrCut]=ExtractVariables(varExtract,param);
    
    fill=[volumefraction(:).targetfill];
    currFill=[volumefraction(:).volumefraction];
    if fillErrCut<2
        fillErrCut=2;
        disp('fillErrCut was too small it was reset to 2')
    end
    newFill=fill;
    if nStep<fillErrStep
        condVals=conds(currFill,fillErrCut);
        newFill=min(newFill,min(min(condVals{1},condVals{2}),1)*0.8+0.2);
        newFill=max(newFill,max(max(condVals{3},condVals{4}),0)*0.8);
        newFill(fill<=0)=-0.1;
        newFill(fill>=1)=1.1;
        for ii=1:numel(volumefraction)
            volumefraction(ii).targetfill=newFill(ii);
        end
    end
    
end


function [dtMax,isConvergingPast,isChangeddtMax]=TestConvergenceRate(resArray,targConv,dtMax,...
        convDistance,isConvergingPast)
    
    resArray=[(1:length(resArray))',log10(resArray)'];
    
    isConverging=true;
    meanRes=mean(resArray(:,2));
    stdRes=std(resArray(:,2));
    minRes=min(resArray(:,2));
    finRes=resArray(end,2);
    startRes=resArray(end,2);
    lastIter=max(resArray(:,1));
    targConv=log10(targConv);
    
    matLin=[resArray(:,1),ones(size(resArray(:,1)))];
    coeffLin=(matLin'*matLin)\matLin'*resArray(:,2);
    linRMS=sqrt(mean((matLin*coeffLin).^2));
    theoretConvIter=round((targConv-coeffLin(2))/coeffLin(1)-lastIter);
    
    if theoretConvIter<0 || theoretConvIter>convDistance
        isConverging=false;
    end
    isChangeddtMax=false;
    if ~isConverging && ~isConvergingPast
        dtMax=dtMax/2;
        isChangeddtMax=true;
    end
    isConvergingPast=isConverging;
end

%% Return Data

function [snakSave]=WriteSnakSave(param,snaxel,dt,snakposition,...
        volumefraction,cellCentredGridSnax,currentConvVelocity,...
        currentConvVolume,movFrame,velcalcinfo,insideContourInfo)
    
    varExtract={'snakData'};
    [snakData]=ExtractVariables(varExtract,param);
    [loopSnaxel]=OrderSurfaceSnaxel(snaxel);
    for ii=1:length(loopSnaxel)
    [xMin(ii),xMax(ii),t(ii),L(ii),A(ii)]=ClosedLoopProperties(loopSnaxel(ii).snaxel.coord);
    end
    
    % Length Snake
    lSnak=sum(L);
    switch snakData
        case 'all'
            [snakSave]=WriteSnakSaveHeavy(snaxel,dt,snakposition,...
                volumefraction,cellCentredGridSnax,currentConvVelocity,...
                currentConvVolume,lSnak,movFrame,velcalcinfo,insideContourInfo);
            
        case 'light'
            [snakSave]=WriteSnakSaveLight(dt,volumefraction,currentConvVelocity,...
                currentConvVolume,lSnak);
            
    end
    
    
    
    
    function [snakSave]=WriteSnakSaveHeavy(snaxel,dt,snakposition,...
            volumefraction,cellCentredGridSnax,currentConvVelocity,...
            currentConvVolume,lSnak,movFrame,velcalcinfo,insideContourInfo)
        
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
        snakSave.lSnak=lSnak;
    end
    
    function [snakSave]=WriteSnakSaveLight(dt,volumefraction,currentConvVelocity,...
            currentConvVolume,lSnak)
        
        snakSave.dt=dt;
        
        snakSave.volumefraction.targetfill=[volumefraction(:).targetfill];
        snakSave.volumefraction.currentfraction=[volumefraction(:).volumefraction];
        snakSave.volumefraction.totVolume=[volumefraction(:).totalvolume];
        
        snakSave.volumefraction.refinedInfo.index=[volumefraction(:).newCellInd];
        snakSave.volumefraction.refinedInfo.fractionvol=[volumefraction(:).splitfraction];
        
        snakSave.currentConvVelocity=currentConvVelocity;
        snakSave.currentConvVolume=currentConvVolume;
        snakSave.lSnak=lSnak;
        
        
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
    
    [edgeLength]=CalculateEgeLengths(unstructured);
    
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
        edgeIndex,loopEdgeIndex,allLoopEdgeIndex,snaxelIndexStart,isInside,edgeLength);
    
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

function [edgeLength]=CalculateEgeLengths(unstructured)
    edgeVertIndex=unstructured.edge.vertexindex';
    vectorEdgeVertIndex=edgeVertIndex(:);
    edgeCoord=unstructured.vertex.coord(FindObjNum([],vectorEdgeVertIndex,unstructured.vertex.index),:);
    
    edgeLength=sqrt(sum((edgeCoord(1:2:end-1,:)-edgeCoord(2:2:end,:)).^2,2));
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
    snaxInd=[snaxel(:).index];
    nextSnaxSub=FindObjNum([],nextSnaxInd,snaxInd);
    precSnaxSub=FindObjNum([],precSnaxInd,snaxInd);
    
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
        edgeIndex,currLoopEdgeIndex,invalidEdgeIndex,snaxelIndexStart,isInside,edgeLength,connectivity)
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
                snaxelIndexStart,initVertexIndex(ii), edgeVertIndex,edgeIndex,edgeLength);
            
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
            generateOrder=FindObjNum([],connecOrder,[snaxelNotOrdered(:).index]);
            if isInside
                generateOrder=generateOrder(end:-1:1);
            end
            [kk2,cellSimVertex{ii},newsnaxel]=GenerateVertexSnaxel...
                (snaxelEdges(generateOrder),kk2,snaxelIndexStart,initVertexIndex(ii),...
                edgeVertIndex,edgeIndex,edgeLength);
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
        snaxelIndexStart,initVertexIndexSingle, edgeVertIndex,edgeIndex,edgeLength)
    
    global snaxInitPos
    
    numSE=length(snaxelEdges); % provides information about snaxel from same vertex
    snaxelEdgesSub=FindObjNum([],snaxelEdges,edgeIndex);
    distRatio=min(edgeLength(snaxelEdgesSub))./edgeLength(snaxelEdgesSub);
    kkLocal=0;
    for jj=1:numSE
        kk=kk+1;
        kkLocal=kkLocal+1;
        snaxIndex=kk+snaxelIndexStart;
        cellSimVertex(jj)=snaxIndex;
        dist=snaxInitPos*distRatio(jj); % Snaxel initialisation, it starts at the vertex
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
    snaxel.orderedge=0; 
end

%% Snaxel Repopulation

function [snaxel,insideContourInfo,nonBreedVert]=...
        SnaxelRepopulate(unstructured,snaxel,insideContourInfo,nonBreedVert...
        ,edgeDat,mergeTopo)
    
    global unstructglobal
    % adds snaxels once a corner has been reached
    [finishedSnakes]=ArrivalCondition(snaxel,edgeDat);
    edgeSnaxel=[snaxel(:).edge];
    edgeVertIndex=unstructured.edge.vertexindex;
    edgeIndex=unstructured.edge.index;
    % This ugly line removes snaxel from the finished list which had
    % arrived at a non breeding vertex
    if ~isempty(nonBreedVert)
        finishedSnakes=finishedSnakes(FindObjNum([],[snaxel(finishedSnakes).tovertex],unique(nonBreedVert))==0);
    end
    %newInsideEdges=[snaxel(finishedSnakes).edge];
    %insideContourInfo(newInsideEdges)=1;
    
    [snaxel,newInsideEdges,delIndex,nonBreedVert]=RepopIterativeBreedingProcess...
        (snaxel,finishedSnakes,edgeVertIndex,edgeIndex,edgeSnaxel,unstructglobal...
        ,nonBreedVert,edgeDat,mergeTopo);
    
    % Removing from repopulation list snaxels that would hit a vertex
    % which has already bred
    [delIndex]=RemoveIdenticalEntries(delIndex(:));
    snaxel=DeleteSnaxel(snaxel,delIndex);
    
    
    [insideContourInfo]=UpdateInsideContourInfo(insideContourInfo,...
        newInsideEdges,snaxel);
    
    
end

function [snaxel,newInsideEdges,delIndex,nonBreedVert]=RepopIterativeBreedingProcess...
        (snaxel,finishedSnakes,edgeVertIndex,edgeIndex,edgeSnaxel,unstructglobal,nonBreedVert,...
        edgeDat,mergeTopo)
    % Iterative Breeding process for the breeding of edges
    
    newInsideEdges=[];
    kk=0;
    delIndex=[];
    rePopInd=[];
    [edgeLength]=CalculateEgeLengths(unstructglobal);
    while ~isempty(finishedSnakes)
        kk=kk+1;
        breedSub=finishedSnakes(1);
        finishedSnakes(1)=[];
        nonBreedVert=[nonBreedVert,snaxel(breedSub).tovertex];
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
            edgeIndex,snaxel(breedSub).edge,snaxel(breedSub).edge,snaxelIndexStart,false,edgeLength,connec);
        
        previousV=snaxel(breedSub).v;
        previousACC=snaxel(breedSub).acc;
        stepL=snaxelRepop(1).d;
%         for ii=1:length(snaxelRepop)
%             snaxelRepop(ii).d=0;
%             snaxelRepop(ii).v=0;%previousV;
%             snaxelRepop(ii).acc=previousACC;
%         end
        rePopInd=[rePopInd,[snaxelRepop(:).index]];
        [snaxel]=AddSnaxel(snaxel,snaxelRepop);
        % Remove breeding snaxels snaxels
        delIndex(kk,1:2)=[snaxel(breedSub).index,indexDoubled]; %#ok<AGROW>
        
        for jj=1:length(snaxelRepop)
            rmvFinSnakes=(edgeSnaxel(finishedSnakes)==snaxelRepop(jj).edge);
            finishedSnakes(rmvFinSnakes)=[];
        end
        
        %finishedSnakes=finishedSnakes(FindObjNum([],[snaxel(finishedSnakes).tovertex],nonBreedVert)==0);
        
    end
    if kk>0
        [snaxel]=TakePostRepopStep(snaxel,rePopInd,stepL,edgeDat,mergeTopo);
    end
end

function [snaxel]=TakePostRepopStep(snaxel,rePopInd,stepL,edgeDat,mergeTopo)
    
    [snaxel]=UpdateSnaxelEdgeOrder(snaxel);
    
    vSnax=[snaxel(:).v];
    
    [snaxel(:).v]=deal(0);
    snaxInd=[snaxel(:).index];
    [maxStepIndiv]=DirectionScaledMaxStep(snaxel,edgeDat);
    repopSub=FindObjNum([],rePopInd,snaxInd);
    rePopInd=rePopInd(repopSub~=0);
    repopSub=repopSub(repopSub~=0);
    
    [snaxel(repopSub).v]=deal(snaxel(repopSub).d);
    snaxelRepopD=[snaxel(repopSub).d];
    [snaxel(repopSub).d]=deal(0);
    
%     for ii=1:numel(repopSub)
%         snaxel(repopSub(ii)).v=(stepL*maxStepIndiv(repopSub(ii)));
%     end
    
    maxDist=MaxTravelDistance(snaxel,mergeTopo);
    
    for ii=1:length(repopSub)
        snaxel(repopSub(ii)).d=max(min([maxDist(repopSub(ii)),snaxelRepopD(ii)]),0);
    end
    for ii=1:length(snaxel)
        snaxel(ii).v=vSnax(ii);
    end
end

function [snaxel]=UpdateSnaxelEdgeOrder(snaxel)
    
    snaxEdge=[snaxel(:).edge];
    snaxOrder=[snaxel(:).orderedge];
    snaxToV=[snaxel(:).tovertex];
    snaxFromV=[snaxel(:).fromvertex];
    
    % This extracts only the edges where there are multiple snaxels on one
    % edge.
    [edgeSnaxList,firstOccurEdge]=unique(snaxEdge,'first');
    snaxEdgeDuplicate=snaxEdge;
    snaxEdgeDuplicate(firstOccurEdge)=[];
    snaxEdgeDuplicate=unique(snaxEdgeDuplicate);
    
    for ii=1:numel(snaxEdgeDuplicate)
        
        posSnaxSub=find(snaxEdgeDuplicate(ii)==snaxEdge);
        snaxOrdAct=snaxOrder(posSnaxSub);
        if any(snaxOrdAct==0)
            unOrdLog=snaxOrdAct==0;
            snaxToSet=snaxToV(posSnaxSub(unOrdLog));
            snaxFromSet=snaxFromV(posSnaxSub(unOrdLog));
            
            newOrders=(snaxFromSet<snaxToSet)+(snaxToSet<snaxFromSet)*(max(snaxOrdAct)+numel(snaxFromSet));
            
            if numel(newOrders)>2
                error('Unexpected number of unclassified snaxels on a single edge')
            end
            if any(newOrders==1)
                snaxOrdAct=snaxOrdAct+1;
            end
            snaxOrdAct(unOrdLog)=newOrders;
            snaxOrder(posSnaxSub)=snaxOrdAct;
        end
    end
    snaxOrder(snaxOrder==0)=1;
    for ii=1:numel(snaxel)
        snaxel(ii).orderedge=snaxOrder(ii);
    end
end

function [snaxel]=TakePostRepopStepOLD(snaxel,rePopInd,stepL,edgeDat)
    vSnax=[snaxel(:).v];
    
    [snaxel(:).v]=deal(0);
    snaxInd=[snaxel(:).index];
    [maxStepIndiv]=DirectionScaledMaxStep(snaxel,edgeDat);
    repopSub=FindObjNum([],rePopInd,snaxInd);
    rePopInd=rePopInd(repopSub~=0);
    repopSub=repopSub(repopSub~=0);
    
    for ii=1:numel(repopSub)
        snaxel(repopSub(ii)).v=(stepL*maxStepIndiv(repopSub(ii)));
    end
    
    maxDist=MaxTravelDistance(snaxel,mergeTopo);
    
    for ii=1:length(repopSub)
        snaxel(repopSub(ii)).d=min([maxDist(repopSub(ii)),stepL]);
    end
    for ii=1:length(snaxel)
        snaxel(ii).v=vSnax(ii);
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
    
    subConnection=FindObjNum([],connection,[snaxel(:).index]);
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

function [snaxel]=AddOrderEdgeField(snaxel)
    % Function to add orderedge to preexisting snaxel restarts.
    
     
    snaxEdge=[snaxel(:).edge];
    snaxD=[snaxel(:).d];
    snaxToV=[snaxel(:).tovertex];
    snaxFromV=[snaxel(:).fromvertex];
    snaxOrder=zeros(size(snaxel));
    % This extracts only the edges where there are multiple snaxels on one
    % edge.
    [edgeSnaxList,firstOccurEdge]=unique(snaxEdge,'first');
    
    
    for ii=1:numel(edgeSnaxList)
        
        actSnaxSub=find(edgeSnaxList(ii)==snaxEdge);
        
        rootVert=min(snaxToV(actSnaxSub(1)),snaxFromV(actSnaxSub(1)));
        
        snaxDAct=snaxD(actSnaxSub)*(snaxFromV(actSnaxSub)==rootVert)+(1-snaxD(actSnaxSub))*(snaxFromV(actSnaxSub)==rootVert);
        [~,orderedSnax]=sort(snaxDAct);
        orderedgeAct=1:numel(snaxDAct);
        snaxOrder(actSnaxSub)=orderedgeAct(orderedSnax);
    end
    
    for ii=1:numel(snaxel)
        snaxel(ii).orderedge=snaxOrder(ii);
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
        
        [contourStruct(ii).vector]=ContourEdgeNormal(snakposition,...
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
%     snaxInd=[snaxel(:).index]';
%     snaxConnects=vertcat(snaxel(:).connectivity);
%     snaxContour1=[snaxInd,snaxConnects(:,1)];
%     snaxContour2=[snaxInd,snaxConnects(:,2)];
%     snaxContour=[snaxContour1;snaxContour2];
%     cellSimilar=FindIdenticalVector(snaxContour);
%     
%     arraySimilar=vertcat(cellSimilar{:});
%     snaxContour=snaxContour(arraySimilar(:,1),:);
%     for ii=length(snaxContour(:,1)):-1:1
%         contourStruct(ii).index1=snaxContour(ii,1);
%         contourStruct(ii).index2=snaxContour(ii,2);
%     end

    [contourStruct(1:numel(snaxel)).index1]=deal(snaxel(:).index);
    [contourStruct(1:numel(snaxel)).index2]=deal(snaxel(:).snaxnext);

end

function [normalVector]=ContourEdgeNormal(snakposition,indices,snaxIndPos)
    % Calculates the normal of the contour
    
    snaxSubPos=FindObjNum([],indices,snaxIndPos);
    
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

%% Velocity Calculation (External Function)

function [snaxel,snakposition,snaxelmodvel,velcalcinfo,forceparam,sensSnax]=VelocityCalculationVolumeFraction...
        (snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    velcalcinfo=[];
    sensSnax=[];
    switch forceparam.velType
        case 'default'
            [snaxel,snakposition,snaxelmodvel,velcalcinfo,sensSnax,forceparam]=...
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
            [snaxel,snakposition,snaxelmodvel,velcalcinfo,sensSnax,forceparam]=...
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

function [finishedSnakesSub]=ArrivalCondition(snaxel,edgeDat)
    % Calculates the arrival condition for repopulation
    
    global arrivalTolerance maxDt snaxInitPos
    v=[snaxel(:).v];
    d=[snaxel(:).d];
    fromvertSnax=[snaxel(:).fromvertex];
    tovertSnax=[snaxel(:).tovertex];
    ordSnax=[snaxel(:).orderedge];
    edgeSnax=[snaxel(:).edge];
    
    %[isClosest]=IsClosestToEndVertex(edgeSnax,fromvertSnax,d);
    [isClosest]=IsClosestToEndVertexOrder(edgeSnax,ordSnax,fromvertSnax,tovertSnax,d);
    [maxStepIndiv]=DirectionScaledMaxStep(snaxel,edgeDat);
    
    isFreeze=logical([snaxel(:).isfreeze]); % finds all unfrozen
    isFwd=v>=0;
    isArriving=d>=(1-arrivalTolerance*maxStepIndiv);
    isArrived=d>=(1-snaxInitPos/100);
    willImpact=abs((1-d)./v)<maxDt;
    finishedSnakesSub=find((isFwd & ~isFreeze & isArriving & willImpact & isClosest)...
        | (isArrived & ~isFreeze & isFwd & isClosest));
end

function [isClosest]=IsClosestToEndVertex(edgeSnax,fromvertSnax,d)
    
    [dGrid1,dGrid2]=meshgrid(d);
    [edgeGrid1,edgeGrid2]=meshgrid(edgeSnax);
    [vertGrid1,vertGrid2]=meshgrid(fromvertSnax);
    logicalMask=(edgeGrid1 == edgeGrid2)-2*((edgeGrid1 == edgeGrid2) & ~(vertGrid1 == vertGrid2));
    logicalMask(logical(eye(numel(d))))=0;
    isClosest=~any((dGrid1-((logicalMask==-1)+logicalMask.*dGrid2))<=0);
end


function [isClosest]=IsClosestToEndVertexOrder(edgeSnax,ordSnax,fromvertSnax,tovertSnax,d)
    
    
    
    [edgeSnaxList,~,posEdgeSnax]=unique(edgeSnax);
    maxOrd=max(ordSnax);
    posEdgeSnax=posEdgeSnax';
    ordExtractMat=nan([numel(edgeSnaxList),maxOrd]); % use nan to ensure min and max not skewed
    indExtractMat=zeros([numel(edgeSnaxList),maxOrd]);
    
    ordExtractMat(sub2ind([numel(edgeSnaxList),maxOrd],posEdgeSnax,ordSnax))=ordSnax;
    indExtractMat(sub2ind([numel(edgeSnaxList),maxOrd],posEdgeSnax,ordSnax))=1:numel(ordSnax);
    
    [~,isMax]=max(ordExtractMat,[],2);
    [~,isMin]=min(ordExtractMat,[],2);
    flagOrderMin=zeros(size(d));
    flagOrderMax=zeros(size(d));
    flagOrderMin(indExtractMat(sub2ind([numel(edgeSnaxList),maxOrd],1:size(indExtractMat,1),isMin')))=-1;
    flagOrderMax(indExtractMat(sub2ind([numel(edgeSnaxList),maxOrd],1:size(indExtractMat,1),isMax')))=1;
    flagLocation=(~xor(d>=0.5,fromvertSnax<tovertSnax))*2-1;
    
    isClosest=(flagOrderMin==flagLocation) | (flagOrderMax==flagLocation);
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

function [snaxel]=InitisaliseSnaxelOrder(snaxel)
   % Initialises the order of all snaxels 
   % Only valid just after Initialisation
   
   snaxFrom=[snaxel(:).fromvertex];
   snaxTo=[snaxel(:).tovertex];
   snaxD=[snaxel(:).d];
   
   snaxOrd=1+xor(((snaxTo-snaxFrom)<0),logical(round(snaxD)));
   for ii=1:numel(snaxel)
       snaxel(ii).orderedge=snaxOrd(ii);
   end
   
end
%% Time Update

function maxDist=MaxTravelDistance(snaxel,mergeTopo)
    % Calculates the maximum distance that can be travelled by a snaxel
    dSnax=[snaxel(:).d];
    vSnax=[snaxel(:).v];
    orderSnax=[snaxel(:).orderedge];
    fromvertSnax=[snaxel(:).fromvertex];
    toVertSnax=[snaxel(:).tovertex];
    edgeSnax=[snaxel(:).edge];
    nSnax=length(edgeSnax);
    if ~mergeTopo
        eps=0;%1e-6; 
    else
        eps=0;
    end
    % here to handle edge cases : avoid 2 points to be next to each other
    % and there order not being well defined.
    maxDist=ones(1,nSnax);
    
    for ii=1:nSnax
        sameEdgeSnax=find(edgeSnax==edgeSnax(ii));
        sameEdgeSnax(sameEdgeSnax==ii)=[];
        
        if vSnax(ii)>=0
            maxDist(ii)=1-dSnax(ii);
        else
            maxDist(ii)=-dSnax(ii);
        end
%         if numel(sameEdgeSnax)==1
%             
%             dSnaxOther=dSnax(sameEdgeSnax);
%             vSnaxOther=vSnax(sameEdgeSnax);
%             sameDir=fromvertSnax(ii)==fromvertSnax(sameEdgeSnax);
%             if sameDir
%                 %disp('Turning Snaxel for impact distance calculation')
%                 dSnaxOther=1-dSnaxOther;
%                 vSnaxOther=-vSnaxOther;
%             end
%             
%             impactDist=(1-dSnax(ii)-dSnaxOther)/...
%                 (vSnax(ii)+vSnaxOther)*vSnax(ii);
%             
%             if vSnax(ii)+vSnaxOther==0
%                 impactDist=1;
%             end
%             
%             if vSnax(ii)>=0 && impactDist>=0
%                 maxDist(ii)=min([maxDist(ii),impactDist-eps]);
%             elseif vSnax(ii)<0 && impactDist<=0
%                 maxDist(ii)=max([maxDist(ii),impactDist+eps]);
%             end
%             
%         else
        if numel(sameEdgeSnax)>0
%             disp('More than 2 snaxels on the same edge')
            
            impactDist=[];
            for jj=1:numel(sameEdgeSnax)
                [impactDist(jj)]=CalcSnaxelRelativeImpact(dSnax(ii),vSnax(ii),fromvertSnax(ii),...
                    dSnax(sameEdgeSnax(jj)),vSnax(sameEdgeSnax(jj)),fromvertSnax(sameEdgeSnax(jj)));
            end
            
            if any(impactDist==0)
                ordersAct=[orderSnax(ii),orderSnax(sameEdgeSnax)];
                [~,orderSortI]=sort(ordersAct);

                orderImpact=sign(toVertSnax(ii)-fromvertSnax(ii))*sign(orderSortI(2:end)-orderSortI(1));
                % + indicates in the way on positive travel
                % - indicates in the way on negative travel
                orderImpact(impactDist~=0)=[];
                if vSnax(ii)>=0 && any(orderImpact>0)
                    maxDist(ii)=0;
                elseif vSnax(ii)<0 && any(orderImpact<0)
                    maxDist(ii)=0;
                end
            end
            if vSnax(ii)>=0 && any(impactDist>0)
                maxDist(ii)=min([maxDist(ii),min(impactDist(impactDist>0))-eps]);
            elseif vSnax(ii)<0 && any(impactDist<0)
                maxDist(ii)=max([maxDist(ii),max(impactDist(impactDist<0))+eps]);
            end
        end
    end
    
end

function [impactDist]=CalcSnaxelRelativeImpact(rootD,rootV,rootVert,otherD,...
        otherV,otherVert)
    
    sameDir=rootVert==otherVert;
    if sameDir
        %disp('Turning Snaxel for impact distance calculation')
        otherD=1-otherD;
        otherV=-otherV;
    end
    vRatio=rootV/(rootV+otherV);
    impactDist=(1-rootD-otherD)*vRatio;
    
    if isnan(impactDist)
        impactDist=0;
    end
    if isinf(impactDist)
        impactDist=sign(impactDist)*1;
    end
end

function [dt,dtSnax2,maxDist]=TimeStepCalculation(snaxel,maxStep,maxDt,dtMin,...
        stepType,vSwitch,edgeDat,mergeTopo)
    % Calculates the timestep to ensure that no snaxel crosses the line without
    % crossing an intersection
    
    vMax=max(abs([snaxel(:).v]));
    
    if vMax<vSwitch
        maxDt=1;
    end
    
    canMove=[snaxel(:).isfreeze]==0;
    maxDist=MaxTravelDistance(snaxel,mergeTopo);
    [maxStepIndiv]=DirectionScaledMaxStep(snaxel,edgeDat);
    maxDist(maxDist>(maxStep*maxStepIndiv))=maxStep*maxStepIndiv(maxDist>(maxStep*maxStepIndiv));
    maxDist=-maxDist;
    maxDist(maxDist>(maxStep*maxStepIndiv))=maxStep*maxStepIndiv(maxDist>(maxStep*maxStepIndiv));
    maxDist=-maxDist;
    
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
    switch stepType
        case 'strict'
            [dt]=TimeStep_Strict(dtStep,dtMax,maxDt,dtSnax,dtMin);
        case 'bounded'
            [dt]=TimeStep_bounded(dtStep,dtMax,maxDt,dtSnax,dtMin);
        case 'mixed'
            [dt]=TimeStep_bounded(dtStep,dtMax,maxDt,dtSnax,dtMin);
            [dtSnax2]=TimeStep_indiv(dtSnax2,dt);
        case 'indiv'
            [dtSnax2]=TimeStep_indiv(dtSnax2,maxDt);
            dt=min(dtSnax2)+0.0001;
    end
    
    if dt<0
        warning('Time going back')
    end
end

function [maxStepIndiv]=DirectionScaledMaxStep(snaxel,edgeDat)
   % scales the maximum step depending on the assymetry of the grid
   edgeOrient=edgeDat.orient;
   edgeIndList=edgeDat.indList;
   refGridRatio=edgeDat.ratio;
   edgeSnax=[snaxel(:).edge];
   maxStepIndiv=ones(size(snaxel));
   maxStepIndiv(logical(edgeOrient(FindObjNum([],edgeSnax,edgeIndList))))=1/refGridRatio;
   maxStepIndiv=maxStepIndiv/max(maxStepIndiv);
end

function [dt]=TimeStep_Strict(dtStep,dtMax,maxDt,dtSnax,dtMin)
    dt=min([dtStep, dtMax, maxDt]);
    if dt==0
        fprintf(' - Breeding Snaxel Ignored - ')
        dt=min([dtSnax(dtSnax~=0),dtMin]);
        if isempty(dt)
            dt=dtMin;
            fprintf(' - Dt forced to DtDefault - ')
        end
    end
end

function [dt]=TimeStep_bounded(dtStep,dtMax,maxDt,dtSnax,dtMin)
    dt=min([dtStep, dtMax, maxDt]);
    if dt<dtMin
        dt=dtMin;
    end
end

function [dtSnax]=TimeStep_indiv(dtSnax,dtMax)
    dtSnax(dtSnax>dtMax)=dtMax;
end

function snaxel=SnaxelDistanceUpdate(snaxel,dt,dtSnax,maxDist,stepType)
    % Updates the current distance of the snaxel
    
    if strcmp(stepType,'indiv') || strcmp(stepType,'mixed')
        snaxel=TimeInaccurateDistanceUpdate(snaxel,dtSnax,maxDist);
    else
        snaxel=TimeAccurateDistanceUpdate(snaxel,dt,maxDist);
    end
end

function snaxel=TimeAccurateDistanceUpdate(snaxel,dt,maxDist)
    count=0;
    
    for ii=1:length(snaxel)
        movDist=snaxel(ii).v*dt;
        
        if abs(movDist)>abs(maxDist(ii))
            count=count+1;
            movDist=maxDist(ii);
        end
        snaxel(ii).d=movDist+snaxel(ii).d;
    end
    if count>0
        fprintf(' - %i Problems with snaxel update',count);
    end
end

function snaxel=TimeInaccurateDistanceUpdate(snaxel,dtSnax,maxDist)
    count=0;
    for ii=1:length(snaxel)
        movDist=snaxel(ii).v*dtSnax(ii);
        if abs(movDist)>abs(maxDist(ii))
            if (abs(movDist)-abs(maxDist(ii))) > 1e-15
                count=count+1;
            end
            movDist=maxDist(ii);
        end
        snaxel(ii).d=movDist+snaxel(ii).d;
    end
    if count>0
        fprintf(' - %i Problems with snaxel update ',count);
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
        canMove=true(size(snaxel)); %[snaxel(:).isfreeze]~=1;
        if ~isempty(snaxel(canMove))
            [delIndex]=FindBadSnaxels(snaxel(canMove),insideContourInfo);
        else
            delIndex=[];
        end
        delSub=FindObjNum([],delIndex,[snaxel(:).index]);
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
                disp('Topology Collapse')
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

function [snaxel]=FreezingFunction(snaxel,borderVertices,edgeDat,mergeTopo)
    % Function which selects snaxels that need to be frozen
    
    if ~exist('mergeTopo','var'),mergeTopo=false;end
    
    if mergeTopo
        edgeFreeze=[];
    else
        [edgeFreeze]=FreezeEdgeContact(snaxel);
    end
    
    [borderFreeze]=FreezeBorderContact(snaxel,borderVertices.strong,edgeDat);
    if isfield(borderVertices,'nonBreedVert')
        %[nonBreedFreeze]=FreezeNonBreedContact(snaxel,borderVertices.nonBreedVert);
        % For this to work proper freeze handling is needed at Velocity
        % level: This means a 2 step process where velocities are based off
        % the actual freeze status after thawing test (ie is it still frozen)
        nonBreedFreeze=[];
        else
        nonBreedFreeze=[];
    end
    if isfield(borderVertices,'nonbreedpers')
%         if ~isempty(borderVertices.nonbreedpers)
            [nonBreedFreeze2]=FreezeNonBreedContact(snaxel,borderVertices.nonbreedpers);
%         else
%             nonBreedFreeze2=[];
%         end
        % For this to work proper freeze handling is needed at Velocity
        % level: This means a 2 step process where velocities are based off
        % the actual freeze status after thawing test (ie is it still frozen)
        nonBreedFreeze=[nonBreedFreeze,nonBreedFreeze2];
    else
        
    end
    
    freezeIndex=[edgeFreeze,borderFreeze];
    
    [freezeIndex]=RemoveIdenticalEntries(freezeIndex);
    [snaxel(:).isfreeze]=deal(false);
    if ~isempty(freezeIndex)
        freezeSub=FindObjNum([],freezeIndex,[snaxel(:).index]);
        for ii=1:length(freezeIndex)
            snaxel(freezeSub(ii)).isfreeze=1;
        end
    end
    if ~isempty(nonBreedFreeze)
        freezeSub=FindObjNum([],nonBreedFreeze,[snaxel(:).index]);
        for ii=1:length(nonBreedFreeze)
            snaxel(freezeSub(ii)).isfreeze=2 ...
            -0.1+round(snaxel(freezeSub(ii)).d)*0.2; 
        % this indicates if the snaxel is blocked at 0 or 1
        % 0 -> No negative speeds -> 1.9 
        % 1 -> No positive speeds -> 2.1
        
        end
    end
end

function [nonBreedFreeze]=FreezeNonBreedContact(snaxel,nonBreedVert)
    
    global arrivalTolerance
    tol=0*arrivalTolerance/1000;
    nonBreedFreeze=[];
    dSnax=[snaxel(:).d];
    vSnax=[snaxel(:).v];
    
    tovertisNonBreed=(FindObjNum([],[snaxel(:).tovertex],unique(nonBreedVert))~=0);
    fromvertisNonBreed=(FindObjNum([],[snaxel(:).fromvertex],unique(nonBreedVert))~=0);
    
    isNonBreedFreeze=(dSnax<=tol & fromvertisNonBreed' & vSnax<0) | ...
        ((1-dSnax)<=tol & tovertisNonBreed' & vSnax>0);
    nonBreedFreeze=[snaxel(isNonBreedFreeze).index];
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
%             if isFreeze(ii)
%                 isImpact(ii)=false;
%             end
            
        end
    end
    freezeIndex=indexSnax(isImpact);
end

function [isImpact]=EdgeImpactCondition(dSnax,vSnax,fromvertSnax,sub1,sub2)
    % edge Impact condition
    
    global snaxInitPos
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
    
    isClose=abs(deltaD)<=snaxInitPos;
    isApproaching=(deltaV)>0;
    isImpact=(isClose && isApproaching) || isContact;
end

function [freezeIndex]=FreezeBorderContact(snaxel,borderVertices,edgeDat)
    % Function which returns the index of Snaxels which have hit the edge of
    % the design space
    
    
    [arrivedSub]=ArrivalCondition(snaxel,edgeDat);
    
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

function [borderVertices]=FindWeakBorderVertex(unstrucReshape,cellCentredGrid,connectstructinfo)
    
    edgeVertIndex=vertcat(unstrucReshape.edge(:).vertexindex);
    edgeIndex=[unstrucReshape.edge(:).index];
    cellIndex=[unstrucReshape.cell(:).index];
    vertexIndex=[unstrucReshape.vertex(:).index];
    edgeCellIndex=vertcat(unstrucReshape.edge(:).cellindex);
    
    [rowPos,colPos]=find(edgeCellIndex==0);
    subCellEdge=FindObjNum([],unique(edgeCellIndex(sub2ind(size(edgeCellIndex),...
        rowPos,abs(colPos-3)))),[connectstructinfo.cell(:).oldCellInd]);
    indFineCell=[connectstructinfo(:).cell(subCellEdge).newCellInd];
    
    subFineCellBord=FindObjNum([],indFineCell,[cellCentredGrid(:).index]);
    vertBord=[cellCentredGrid(unique(subFineCellBord)).vertex];
    borderVertices=[vertBord(:).index];
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
    
    workSub=FindObjNum([],workPair,[snaxel(:).index]);
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

function [freezeVertex,borderEdgInd]=IdentifyProfileEdgeVertices(refinedGrid)
    
    borderEdges=[refinedGrid.edge(:).boundaryis0] & [refinedGrid.edge(:).boundaryis1];
    borderEdgInd=[refinedGrid.edge(borderEdges).index];
    borderVertices=sort([refinedGrid.edge(borderEdges).vertexindex]);
    freezeVertex=RemoveIdenticalEntries(borderVertices);
    
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
                
                % coord=oldGrid.cell(oldCellIndUnstructSub(ii)).coord;
                % frac=volumefraction(oldCellIndVolFracSub(ii)).volumefraction...
                %     -volumefraction(oldCellIndVolFracSub(ii)).targetfill;
                %frac=volumefraction(oldCellIndVolFracSub(ii)).oldCellInd;
                % PlotVolFrac(figh,axh,coord,frac)
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

function [figh]=CheckResultsLight(unstructured,snakposition,snaxel,figh)
    global nDim domainBounds
    
    if nDim==2
        if nargin==3
            figh=figure;
        else
            figure(figh)
        end
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
    if numel(snaxel(1).v)==1
        for ii=1:length(snakposition)
            X(ii)=snakposition(ii).coord(1);
            Y(ii)=snakposition(ii).coord(2);
            U(ii)=snakposition(ii).vector(1)*snaxel(ii).v/40;
            V(ii)=snakposition(ii).vector(2)*snaxel(ii).v/40;
        end
    else
        for ii=1:length(snakposition)
            X(ii)=snakposition(ii).coord(1);
            Y(ii)=snakposition(ii).coord(2);
            U(ii)=snaxel(ii).v(1)/1000;
            V(ii)=snaxel(ii).v(2)/1000;
        end
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
        U(ii)=0*snakposition(ii).vector(1)/80;
        V(ii)=0*snakposition(ii).vector(2)/80;
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
        plot(X,Y,'o--')
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
        dChange{ii}=sensSnax(:,ii)/max(abs(sensSnax(:,ii)))*1;
        %dChange{ii}=-sensSnax(:,ii)/100;
        dAct=dCurr'+dChange{ii};
        
        for jj=1:length(snaxel)
            snaxCopy(jj).d=dAct(jj);
        end
        [snakposition2]=PositionSnakes(snaxCopy,unstructglobal);
        
        
        coord2=vertcat(snakposition2(:).coord);
        Delta{ii}=coord2-coord1;
        [testPos2]=CreateComparisonMatrix(coord2);
        [newOrd]=CompareTestpos(testPos1,testPos2,snaxOrd,dir);
        %newOrd=snaxOrd;
        figure
        plot(coord1(snaxOrd,1),coord1(snaxOrd,2),'+-',coord2(newOrd,1),coord2(newOrd,2),'o-')
        hold on
        for jj=1:length(newOrd)
            plot([coord1(newOrd(jj),1),coord2(newOrd(jj),1)],[coord1(newOrd(jj),2),coord2(newOrd(jj),2)],'k--')
        end
        title(['mode ',int2str(ii)])
        
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
    newOrd=oldOrd(newOrd);
end

%% Various

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
