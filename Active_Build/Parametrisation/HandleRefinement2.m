%% Refined Optimisation steps

function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        HandleRefinement(paramoptim,iterstruct,outinfo,oldBase,gridrefined,...
        connectstructinfo,refStep,nIter,firstValidIter)
    % This must only recieve a portion of iterstruct
    
    oldGrid.base=oldBase;
    oldGrid.refined=gridrefined;
    oldGrid.connec=connectstructinfo;
    oldGrid.cellrefined=CellCentreGridInformation(gridrefined);
    
    % grid refinement
    [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        InitialiseRefinement(paramoptim,iterstruct,outinfo,oldGrid,refStep);
    
   
    [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,nIter,firstValidIter);
    
end


% 

function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        InitialiseRefinement(paramoptim,iterstruct,outinfoOld,oldGrid,refStep)
    
    procStr='REFINE OPTIMISATION PROCESS';
    [tStart]=PrintStart(procStr,1);
    
    % Initialise Optimisation
    % Get Parametrisation parameters
    varNames={'optimCase'};
    optimCase=ExtractVariables(varNames,paramoptim);
    paramoptim=SetVariables(varNames,{[optimCase,'_',int2str(refStep)]},paramoptim);
    varNames={'boundstr','corneractive','defaultCorner'};
    [boundstr,corneractive,defaultCorner]=ExtractVariables(varNames,paramoptim.parametrisation);
    
    paramoptim.general.desvarconnec=[];
    
    [outinfo]=OptimisationOutput('init',paramoptim);
    warning ('[~,paramoptim]=ConstraintMethod(''init'',paramoptim,[]); Not supported');

    % Refine Grid
    varNames={'refineOptim'};
    refineCellLvl=ExtractVariables(varNames,paramoptim);
    refparamsnake=SetVariables({'refineGrid'},{refineCellLvl(refStep,:)},...
        paramoptim.parametrisation);
    
    [~,baseGrid,gridrefined,connectstructinfo,~,~]...
        =GridInitAndRefine(refparamsnake,oldGrid.base);
    
    newgrid.base=baseGrid;
    newgrid.refined=gridrefined;
    newgrid.connec=connectstructinfo;
    newgrid.cellrefined=CellCentreGridInformation(gridrefined);
    
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation,gridrefined);
    
    % Update some parameters
    
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [paramoptim.general.desvarconnec]=...
        ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
    % Update Fill Information to match snake
    varNames={'refineGrid'};
    refineGrid=ExtractVariables(varNames,paramoptim.parametrisation);
    if numel(refineGrid)==1; refineGrid=ones(1,2)*refineGrid;end
    [gridmatch,~]=GridMatching(oldGrid,newgrid,refineGrid,refineCellLvl(refStep,:));
    
        % Fill follows the order of activeCells
    
    [profloops]=ExtractVolInfo(outinfoOld.rootDir);
    
    [profloops,transformstruct]=ExtractVolumeFraction(gridmatch,profloops);
    
    % Generate new snake restarts.
    [baseGrid,gridrefined]=ReFracGrids(baseGrid,gridrefined,...
        connectstructinfo,profloops(1).newFracs);
    
    if ~corneractive
        [baseGrid,gridrefined]=ReduceCornerFrac(baseGrid,gridrefined,...
        connectstructinfo,defaultCorner);
    end
    
    [gridrefined]=EdgePropertiesReshape(gridrefined);
    [loop]=GenerateSnakStartLoop(gridrefined,boundstr);
    
    iterstruct=RewriteHistory(iterstruct,profloops,baseGrid);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,...
        paramoptim.spline,outinfo,0,0,0);
    
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    
    [~]=PrintEnd(procStr,1,tStart);
end

function [profloops]=ExtractVolInfo(optimRootDir)
    
    [iterationPath,iterName]=FindDir(optimRootDir,'iteration',true);
    iterNum=regexp(iterName,'iteration_','split');
    
    profloops=repmat(struct('iter',[],'prof',[],'refinevolfrac',[]),[1 0]);
    
    for ii=1:length(iterationPath)
        
        [profloops]=[profloops,FindProfileLoops(iterationPath{ii},...
            str2double(iterNum{ii}{2}))];
    end
    
end

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
end

function [profloops]=FindProfileLoops(rootDir,iterNum)
    
    [profilePath,profileName]=FindDir(rootDir,'profile',true);
    profNum=regexp(profileName,'profile_','split');
    %profNum=profNum(:,2);
    kk=1;
    for ii=1:length(profilePath)
        
        [restartPath,restartName]=FindDir(profilePath{ii},'restart',false);
        load(restartPath{1},'snakSave')
        
        profloops(kk).iter=iterNum;
        profloops(kk).prof=str2num(profNum{ii}{2});
        profloops(kk).refinevolfrac=snakSave(end).volumefraction.refinedInfo; %#ok<COLND>
        kk=kk+1;
    end
    
    
    
end

function [profloops,transformstruct]=ExtractVolumeFraction(gridmatch,profloops)
    
    [transformstruct,~]=BuildMatrix(gridmatch);
    
    [profloops]=ConvertProfToFill(profloops,transformstruct);
end

function [transformstruct,coeffMat]=BuildMatrix(gridmatch)
    
    newGridIndsMulti=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {gridmatch.matchstruct(:).oldGridInd},...
        {gridmatch.matchstruct(:).newGridInd},'UniformOutput',false));
    oldGridInds=[gridmatch.matchstruct(:).oldGridInd];
    oldGridCoeff=[gridmatch.matchstruct(:).coeff];
    newGridInds=[gridmatch.matchstruct(:).newGridInd];
    
    oldGridUniq=RemoveIdenticalEntries(sort(oldGridInds));
    oldGridSub=FindObjNum([],oldGridInds,oldGridUniq);
    newGridMultiSub=FindObjNum([],newGridIndsMulti,newGridInds);
    
    coeffMat=zeros([numel(newGridInds),numel(oldGridUniq)]);
    coeffMat(sub2ind(size(coeffMat),newGridMultiSub,oldGridSub))=oldGridCoeff;
    transformstruct.coeff=coeffMat;
    transformstruct.indNew=newGridInds;
    transformstruct.indOld=oldGridUniq;
    transformstruct.volumeNew=[gridmatch.matchstruct(:).newvolume]';
end

function [profloops]=ConvertProfToFill(profloops,transformstruct)
    for ii=1:length(profloops)
        volSubs=FindObjNum([],transformstruct.indOld,profloops(ii).refinevolfrac.index);
        profloops(ii).newFracs=(transformstruct.coeff*...
            profloops(ii).refinevolfrac.fractionvol(volSubs)')./ transformstruct.volumeNew;
    end
end


function [newGrid,newRefGrid]=ReFracGrids(baseGrid,refinedGrid,...
        connectstructinfo,newFracs)
    
    activeCell=true(size(baseGrid.cell));
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    if numel(newFracs)~=numel(activeCellSub)
        error('Fill and Active Set do not match in size')
    end
    if sum(abs(newFracs))==0
        newFracs(round(numel(newFracs)/2))=1e-3;
    end
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFracs(ii);
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFracs(ii));
    end
    
end


function [newGrid,newRefGrid]=ReduceCornerFrac(baseGrid,refinedGrid,...
        connectstructinfo,cornerFill)
    
    activeCell=~logical([baseGrid.cell(:).isactive]);
    activeFill=logical([baseGrid.cell(:).fill]~=0);
    activeCell=activeCell & activeFill;
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=cornerFill;
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(cornerFill);
    end
    
end

function [loop]=GenerateSnakStartLoop(gridrefined2,boundstr)
    
    isEdge=[gridrefined2.edge(:).(boundstr{1})];
    cond=boundstr{3};
    [loop]=OrderSurfaceVertexReshape(gridrefined2,isEdge,cond);
    
end

function [iterstruct]=RewriteHistory(iterstruct,profloops,baseGrid)
    
    actFill=logical([baseGrid.cell(:).isactive]);
    
    for ii=find([profloops(:).iter]~=0)
        
        iterstruct(profloops(ii).iter).population(profloops(ii).prof).fill=...
            profloops(ii).newFracs(actFill)';
        
    end
    
    % rewrite optimdat
    iterstruct(1).population(1).optimdat.value=0;
    iterstruct(1).population(1).optimdat.var=1;
    for ii=3:2:length(iterstruct)
        iterstruct(ii).population(1).optimdat.value=iterstruct(ii).population(1).fill...
            -iterstruct(ii-1).population(1).fill;
        
        iterstruct(ii).population(1).optimdat.var=1:numel(iterstruct(ii-1).population(1).fill);
    end
    for ii=1:length(iterstruct)
        for jj=2:length(iterstruct(ii).population)
            iterstruct(ii).population(jj).optimdat.value=...
                iterstruct(ii).population(jj).fill-iterstruct(ii).population(1).fill;
            iterstruct(ii).population(jj).optimdat.var=find(...
                iterstruct(ii).population(jj).optimdat.value~=0);
            iterstruct(ii).population(jj).optimdat.value=...
                iterstruct(ii).population(jj).optimdat.value(...
                iterstruct(ii).population(jj).optimdat.var);
        end
    end
end














%% From executeoptimisation
%{
function [unstructured,unstructReshape,gridrefined,connectstructinfo,...
        unstructuredrefined,loop]=GridInitAndRefine(param,unstructReshape)
    % Executes the Grid Initialisation process
    procStr='INITIAL GRID OPERATIONS';
    [tStart]=PrintStart(procStr,2);
    
    if nargin<2
        [unstructured,~,unstructReshape]=GridInitialisationV2(param);
    else
        [unstructured]=ModifReshape(unstructReshape);
    end
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    
    
    [~]=PrintEnd(procStr,2,tStart);
    
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

function [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid)
    
    varExtract={'symType'};
    [symType]=ExtractVariables(varExtract,paramoptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    nDesVar=sum([baseGrid.cell(:).isactive]);
    paramoptim.general.nDesVar=nDesVar;
    
    paramoptim.general.symDesVarList...
        =BuildSymmetryLists(symType,cellLevels,corneractive);
    
    paramoptim.general.notDesInd...
        =BuildExclusionList(paramoptim.general.symDesVarList);
    
    [paramoptim]=CheckiterGap(paramoptim);
    
    
    % Set resampleSnak to match iin both with precedence in the
    % optimisation option
    paramoptim.parametrisation=SetVariables({'resampleSnak'},...
        {ExtractVariables({'resampleSnak'},paramoptim)},paramoptim.parametrisation);
    
    
end

function [paramoptim]=CheckiterGap(paramoptim)
    varExtract={'optimMethod'};
    [optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    switch optimMethod
        case 'conjgradls'
            paramoptim.general.iterGap=2;
        case 'conjgrad'
            paramoptim.general.iterGap=2;
        otherwise
            
            paramoptim.general.iterGap=1;
    end
end

function [symDesVarList]=BuildSymmetryLists(symType,cellLevels,corneractive)
    
    switch symType
        case 'none'
            symDesVarList=zeros([2,0]);
        case 'horz'
            warning('Symmetry assignement not robust when taken with the rest of the program')
            nRows=cellLevels(2);
            rowMatch=zeros([2,floor(nRows/2)]);
            for ii=1:floor(nRows/2)
                rowMatch(:,ii)=[ii;(nRows-ii+1)];
            end
            
            % Section to deal with inactive corners
            rowLength=ones([1,nRows])*cellLevels(1);
            if ~corneractive
                rowLength(1)=rowLength(1)-2;
                rowLength(end)=rowLength(end)-2;
            end
            indStart=cumsum([0,rowLength]);
            % Section to build index matching lists
            cellMatch=zeros([2,indStart(floor(nRows/2)+1)]);
            for ii=1:floor(nRows/2)
                for jj=1:rowLength(ii)
                    rootInd=jj+indStart(ii);
                    mirrorInd=jj+indStart(cellLevels(2)-(ii-1));
                    
                    cellMatch(:,rootInd)=[rootInd,mirrorInd];
                end
            end
            
            symDesVarList=cellMatch;
            
        case 'vert'
            error('Not coded yet')
            
    end
    
end

function notDesInd=BuildExclusionList(symDesVarList)
    
    notDesInd=[];
    notDesInd=[notDesInd,symDesVarList(2,:)];
    
    
    notDesInd=sort(RemoveIdenticalEntries(notDesInd));
    
end

function [desvarconnec,cellGrid,vertexGrid]=ExtractVolumeCellConnectivity(unstructReshape)
    % This relies on there being no empty indices for speed
    
    [cellGrid]=CellCentredGrid(unstructReshape);
    [vertexGrid]=VertexCentredGrid(unstructReshape);
    
    nCells=length(unstructReshape.cell);
    nEdges=length(unstructReshape.edge);
    nVerts=length(vertexGrid);
    cellList=[unstructReshape.cell(:).index];
    
    desvarconnec=struct('index',[],'neighbours',[],'corners',[]);
    desvarconnec=repmat(desvarconnec,[1,nCells]);
    for ii=1:nCells
        desvarconnec(ii).index=unstructReshape.cell(ii).index;
    end
    for ii=1:nEdges % Find edge connected cells
        cellInd=unstructReshape.edge(ii).cellindex;
        cellInd(cellInd==0)=[];
        cellSub=FindObjNum([],cellInd,cellList);
        for jj=1:length(cellInd)
            desvarconnec(cellSub(jj)).neighbours=...
                [desvarconnec(cellSub(jj)).neighbours,cellInd([1:jj-1,jj+1:end])];
        end
    end
    for ii=1:nVerts % Find vertex connected cells
        cellInd=[vertexGrid(ii).cell(:).index];
        cellInd(cellInd==0)=[];
        cellSub=FindObjNum([],cellInd,cellList);
        for jj=1:length(cellInd)
            desvarconnec(cellSub(jj)).corners=...
                [desvarconnec(cellSub(jj)).corners,cellInd([1:jj-1,jj+1:end])];
        end
    end
    for ii=1:nCells % Remove duplicates sort and build corners
        desvarconnec(ii).neighbours=unique(desvarconnec(ii).neighbours);
        desvarconnec(ii).corners=unique(desvarconnec(ii).corners);
        [~,~,IB] = intersect(desvarconnec(ii).neighbours,desvarconnec(ii).corners);
        desvarconnec(ii).corners(IB)=[];
    end
    
    
end

function [desvarconnec]=ExtractDesignVariableConnectivity(baseGrid,desvarconnec)
    
    actCell=[baseGrid.cell(:).isactive];
    indCell=[baseGrid.cell(:).index];
    desCellInd=cumsum(actCell).*actCell;
    desIndC{length(actCell)}=[];
    
    actSub=find(actCell);
    
    for ii=actSub
        desIndC{ii}=desCellInd(ii);
    end
    for ii=actSub
        desvarconnec(ii).index=desIndC{ii};
        desvarconnec(ii).neighbours=[desIndC{FindObjNum([],desvarconnec(ii).neighbours,indCell)}];
        desvarconnec(ii).corners=[desIndC{FindObjNum([],desvarconnec(ii).corners,indCell)}];
    end
    desvarconnec=desvarconnec(actSub);
end
%% Print to screen functions

function [tStart]=PrintStart(procStr,lvl)
    
    procStart=[procStr,' start'];
    tStart=now;
    switch lvl
        case 0
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp('********************************************************************************************')
            disp('--------------------------------------------------------------------------------------------')
            disp(procStart)
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
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
        case 0
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp(procStart)
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('********************************************************************************************')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
        case 1
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(procStart)
            disp('--------------------------------------------------------------------------------------------')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
            
        case 2
            
            
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(procStart)
            disp('-----------------------')
            disp('  ')
        case 3
            disp('----------')
            
    end
    
end
%}