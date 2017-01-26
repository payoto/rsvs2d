%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%           - TecPlot Outputs - 
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function []=TecplotOutput(entryPoint,varargin)
    
    switch entryPoint
        case 'snakes'
            TecplotOutput_snakes(varargin{:});
        case 'optiminit'
            TecplotOutput_OptimInit(varargin{:});
        case 'optim'
            TecplotOutput_Optim(varargin{:});
    end
        
    
end

function []=TecplotOutput_snakes(FID,baseGrid,fineGrid,snakSave,connectstructinfo)
    
    [baseCellCentredGrid]=CellCentredGridInformationReduced(baseGrid);
    
    varsCell={};%{['VARIABLES = "X" ,"Y", "U" ,"V", "MAG" ,"TARGFILL" ,"VOLFRAC", "DIFF"']};
    [cellBaseGrid]=GridStructToCellOut(baseGrid,1);
    [cellFineGrid]=GridStructToCellOut(fineGrid,2);
    [cellConv]=ConvergenceStructToCellOut(snakSave,5:8);
    time=0;
    for ii=1:length(snakSave)
        [cellSnax(ii).cells]=SnaxelToCellOut(snakSave(ii).snaxel,snakSave(ii).snakposition,3,time);
        if ii==1
            connShare=[];
            varshare=[];
        else
            connShare=4+length(snakSave);
            varshare.zone=connShare;
            varshare.vars=[1,2];
        end
        [cellVol(ii).cellMesh]=CellCentredMeshDataExtraction(baseGrid,snakSave(ii).volumefraction,4,time,connShare,varshare);
        time=time+snakSave(ii).dt;
        %time=time+1;
    end
    
    [cellMesh]=CellCentredMeshDataExtraction(baseGrid,snakSave(1).volumefraction,[],[],[],[]);
    cellWrite=[varsCell,cellBaseGrid,cellFineGrid,cellMesh,[cellSnax(:).cells],[cellVol(:).cellMesh],cellConv];
    WriteToFile(cellWrite,FID);
    fclose(FID);
end

function []=TecplotOutput_OptimInit(FID,baseGrid,fineGrid,cellCentredGrid,...
        snaxel,snakposition)
    
    [cellBaseGrid]=GridStructToCellOut(baseGrid,1);
    [cellFineGrid]=GridStructToCellOut(fineGrid,2);
    
    [cellMesh]=CellCentredMeshDataExtraction_Optim(baseGrid,cellCentredGrid,[],[],[],[]);
    [cellSnax]=SnaxelToCellOut(snaxel,snakposition,3,0);   
    [cellVol]=CellCentredMeshDataExtraction_Optim(baseGrid,cellCentredGrid,4,...
        0,[],[]);
    
    cellWrite=[cellBaseGrid,cellFineGrid,cellMesh,cellSnax,cellVol];
    
    WriteToFile(cellWrite,FID);
    fclose(FID);
end

function []=TecplotOutput_Optim(FID,baseGrid,cellCentredGrid,snaxel,snakposition,time)
    
    [cellSnax]=SnaxelToCellOut(snaxel,snakposition,3,time);
        connShare=5;
        varshare.zone=connShare;
        varshare.vars=[1,2];
   
    [cellVol]=CellCentredMeshDataExtraction_Optim(baseGrid,cellCentredGrid,4,...
        time,connShare,varshare);
    
    cellWrite=[cellSnax,cellVol];
    
    WriteToFile(cellWrite,FID);
    fclose(FID);
end

function []=TecplotOutput_Final()
    
    % Will need to be built
    [cellConv]=ConvergenceStructToCellOut(snakSave,5:8);
    
end

function [cellMesh]=CellCentredMeshDataExtraction_Optim(gridStruct,cellCentredGrid,strandID,time,connShare,varshare)
    
%     cellcentredstruct=repmat(struct('targetfill',[],'volumefraction',[]),...
%         1,length(cellCentredGrid));
    if numel(cellCentredGrid)==1
        for ii=1:length(cellCentredGrid.targetfill)
                cellcentredstruct(ii).targetfill=cellCentredGrid.targetfill(ii);
                cellcentredstruct(ii).volumefraction=cellCentredGrid.targetfill(ii);
        end
    else
        for ii=1:length(cellCentredGrid)
                cellcentredstruct(ii).targetfill=cellCentredGrid(ii).targetfill;
                cellcentredstruct(ii).volumefraction=cellCentredGrid(ii).volumefraction;
        end
    end
   [cellMesh]=CellCentredMeshDataExtraction(gridStruct,cellcentredstruct,strandID,time,connShare,varshare);
end
%% Open File

function []=WriteToFile(cellData,FID)
    % writes the data to the file
    
    
    for ii=1:length(cellData)
        for jj=1:length(cellData{ii}(:,1))
            fprintf(FID,[cellData{ii}(jj,:),'\n']);
        end
    end
    
end

function [FID]=OpenValidFile(optionalSubFolder,datType)
    % Creates a file in the current directory to write data to.
    
    datestr(now)
    [targetDir,fileName]=OutputDirectory(datType,optionalSubFolder);
    
    FID=fopen([targetDir,'\',fileName],'w+');
    
end

function [targetDir,fileName]=OutputDirectory(datType,optionalSubFolder)
    % Creates output directory and extracts valid file name from date and
    % time
    if ~exist('optionalSubFolder','var'),optionalSubFolder='';end
    
    fileName=['TecPlot360_',datestr(now,30),'_',datType,'.plt'];
    stampedSubFolders=['DataArchive_',datestr(now,'yyyy_mm')];
    targetDir=[cd,'\..\results\',optionalSubFolder,'\',stampedSubFolders,'\'];
   
    mkdir(targetDir)
    
end

%% Tecplot Headers for various zone types

function cellHeader=FELINESEGHeader(numNodes,numElm,strandID,time)
    
%     if ~exist('numNodes','var');numNodes=1;end
%     if ~exist('numElm','var');numElm=1;end
    
    
    nodesStr=['NODES=',int2str(numNodes),' '];
    elmStr=['ELEMENTS=',int2str(numElm)];
    idStr=['STRANDID=',int2str(strandID)];
    
    kk=1;
    
    cellHeader{kk}=['VARIABLES = "X" ,"Y" ,"U" ,"V", "MAG"']; kk=kk+1;
    cellHeader{kk}=['ZONE']; kk=kk+1;
    cellHeader{kk}=['VARLOCATION=([1-5]=NODAL)']; kk=kk+1;
    cellHeader{kk}=nodesStr; kk=kk+1;
    cellHeader{kk}=elmStr; kk=kk+1;
    if exist('time','var');
        cellHeader{kk}=idStr; kk=kk+1;
        timeStr=['SOLUTIONTIME=',num2str(time,'%.24f'),' '];
        cellHeader{kk}=timeStr; kk=kk+1;
    end
    cellHeader{kk}='DATAPACKING=POINT'; kk=kk+1;
    cellHeader{kk}='ZONETYPE=FELINESEG'; kk=kk+1;
    
    
    
end

function cellHeader=CellCentredDataHeader(numNodes,numElm,nFaces,vararg)
    % vararg is a structure that can be used to specify arguments:
    % Possible vararg fields
    % strandID, time, connecShareZone, varsharezone
    
    varFields=fieldnames(vararg);
    for ii=1:length(varFields)
        eval([varFields{ii},'=vararg.(varFields{ii});'])
    end
    
    
    nodesStr=['NODES=',int2str(numNodes),' '];
    elmStr=['ELEMENTS=',int2str(numElm)];
    faceStr=['FACES=',int2str(nFaces)];

    kk=1;
    
    cellHeader{kk}=['VARIABLES = "X" ,"Y" ,"TARGFILL" ,"VOLFRAC", "DIFF"']; kk=kk+1;
    cellHeader{kk}=['ZONE']; kk=kk+1;
    cellHeader{kk}=['VARLOCATION=([1-2]=NODAL ,[3-5]=CELLCENTERED)']; kk=kk+1;
    cellHeader{kk}=nodesStr; kk=kk+1;
    cellHeader{kk}=elmStr; kk=kk+1;
    cellHeader{kk}=faceStr; kk=kk+1;
    cellHeader{kk}=['NUMCONNECTEDBOUNDARYFACES=0']; kk=kk+1;
    cellHeader{kk}=['TOTALNUMBOUNDARYCONNECTIONS=0']; kk=kk+1;
    if ~isempty(time);
        idStr=['STRANDID=',int2str(strandID)];
        cellHeader{kk}=idStr; kk=kk+1;
        timeStr=['SOLUTIONTIME=',num2str(time,'%.24f'),' '];
        cellHeader{kk}=timeStr; kk=kk+1;
    end
    if ~isempty(connecShareZone);
        connecZoneStr=['CONNECTIVITYSHAREZONE=',int2str(connecShareZone),' '];
        cellHeader{kk}=connecZoneStr; kk=kk+1;
    end
    if ~isempty(varsharezone);
        cellHeader{kk}=VariableShareZoneString(varsharezone); kk=kk+1;
    end
    cellHeader{kk}='DATAPACKING=BLOCK'; kk=kk+1;
    cellHeader{kk}='ZONETYPE=FEPOLYGON'; kk=kk+1;
    
end

function [str]=VariableShareZoneString(varsharezone)
    % Creates teh string for Variable sharing depending on the structure
    % varsharezone with fields .vars and .zone
    
    str='VARSHARELIST=(';
    for ii=1:length(varsharezone)
        str=[str,'['];
        for jj=1:length(varsharezone(ii).vars)
            str=[str,int2str(varsharezone(ii).vars(jj))];
            if jj~=length(varsharezone(ii).vars)
                str=[str,','];
            end
        end
        str=[str,']=',int2str(varsharezone(ii).zone)];
        if ii~=length(varsharezone)
            str=[str,','];
        end
    end
    str=[str,')'];
end


%% Extract Data for cell formatting

function [cellMesh]=GridStructToCellOut(gridStruct,strandID)
    
    coordDat=vertcat(gridStruct.vertex(:).coord);
    vertIndex=[gridStruct.vertex(:).index];
    connDat=vertcat(gridStruct.edge(:).vertexindex);
    vectorDat=zeros(size(coordDat));
    velDat=vectorDat(:,1);
    vectorDat=[vectorDat,velDat];
    [cellMesh]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID);
    
end

function [cellMesh]=ConvergenceStructToCellOut(snakSave,strandID)
    
    iter=1:length(snakSave);
    convVolume=log10([snakSave(:).currentConvVolume]);
    convVolume(isinf(convVolume(convVolume<0)))=-1000;
    convVelocity=log10([snakSave(:).currentConvVelocity]);
    
    convVelocity(isinf(convVelocity(convVelocity<0)))=-1000;
    padding=zeros(size(iter));
    
    vertIndex=iter;
    vectorDat=[convVolume',convVelocity',padding'];
    connDat=[[1:length(snakSave)-1]',[2:length(snakSave)]'];
    
    coordDat=[iter',convVolume'];
    [cellMesh1]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID(1));
    time=0;
    for ii=1:length(vertIndex)
        [cellMesh3(ii).cells]=CellEdgeMesh(coordDat([ii],:),1,[1 1],vectorDat([ii],:),strandID(3),time);
        time=time+snakSave(ii).dt;
    end
    
    coordDat=[iter',convVelocity'];
    [cellMesh2]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID(2));
    time=0;
    for ii=1:length(vertIndex)
        [cellMesh4(ii).cells]=CellEdgeMesh(coordDat([ii],:),1,[1 1],vectorDat([ii],:),strandID(4),time);
        time=time+snakSave(ii).dt;
    end
    
    
    cellMesh=[cellMesh1,cellMesh2,[cellMesh3(:).cells],[cellMesh4(:).cells]];
end

function [cellMesh]=SnaxelToCellOut(snaxel,snakposition,strandID,time)
    
    
    coordDat=vertcat(snakposition(:).coord);
    vectorDat=vertcat(snakposition(:).vector);
    velDat=[snaxel(:).v]';
    velDat=[[snaxel(:).isfreeze]']*2+1;
    velDat=[[snaxel(:).orderedge]'];
    vectorDat(:,1)=vectorDat(:,1).*velDat;
    vectorDat(:,2)=vectorDat(:,2).*velDat;
    vectorDat=[vectorDat,velDat];
    vertIndex=[snakposition(:).index];
    connDat=[[snaxel(:).index]',[snaxel(:).snaxnext]'];
    
    [cellMesh]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID,time);
    
end

function [cellMesh]=CellCentredMeshDataExtraction(gridStruct,cellcentredstruct,strandID,time,connShare,varshare)
    % Saves the mesh as a cell-centred polygonal mesh
    
    coordDatHoriz=vertcat(gridStruct.vertex(:).coord)';
    
    vertIndex=[gridStruct.vertex(:).index];
    connNodeDat=vertcat(gridStruct.edge(:).vertexindex);
    connNodeDat(:,1)=FindObjNum([],connNodeDat(:,1),vertIndex);
    connNodeDat(:,2)=FindObjNum([],connNodeDat(:,2),vertIndex);
    
    cellIndex=[gridStruct.cell(:).index];
    connElmDat=vertcat(gridStruct.edge(:).cellindex);
    connElmDat(:,1)=FindObjNum([],connElmDat(:,1),cellIndex);
    connElmDat(:,2)=FindObjNum([],connElmDat(:,2),cellIndex);
    connElmDat=connElmDat';
    
    cCTargFill=[cellcentredstruct(:).targetfill];
    cCVolFrac=[cellcentredstruct(:).volumefraction];
    cCDiff=abs(cCTargFill-cCVolFrac);
    
    cellCentredData=[cCTargFill;cCVolFrac;cCDiff];
    
    vararg.strandID=strandID;
    vararg.time=time;
    vararg.connecShareZone=[];
    if connShare~=0
        vararg.connecShareZone=connShare;
    end
    vararg.varsharezone=[];
    if ~isempty(varshare)
        vararg.varsharezone=varshare;
    end
    
    [cellMesh]=CellCellCentredMesh(coordDatHoriz,cellCentredData,connNodeDat,connElmDat,vararg);
end

%% Format data into cell arrays of strings Line by line

function cellLoops=CellLoop(loopout)
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
                    num2str(loopout.surf(ii).coord(jj,ll),24),'  '];
            end
            kk=kk+1;
        end
        
    end
    
end

function [cellMesh]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID,time)
    
    
    [numNode,nDim]=size(coordDat);
    [numElm,nConn]=size(connDat);
    % Zone Header
    if ~exist('time','var')
        cellHeader=FELINESEGHeader(numNode,numElm,strandID);
    else
        cellHeader=FELINESEGHeader(numNode,numElm,strandID,time);
    end
    % Coordinates
    for ii=numNode:-1:1
        cellNodeCoord{ii}='';
        for jj=1:nDim
            cellNodeCoord{ii}=[cellNodeCoord{ii},...
                num2str(coordDat(ii,jj),24),'  '];
        end
        for jj=1:nDim+1
            cellNodeCoord{ii}=[cellNodeCoord{ii},...
                num2str(vectorDat(ii,jj),24),'  '];
        end
    end
    % Connectivity
    for ii=numElm:-1:1
        cellConnCoord{ii}='';
        for jj=1:nConn
            cellConnCoord{ii}=[cellConnCoord{ii},...
                int2str(FindObjNum([],connDat(ii,jj),vertIndex)),'  '];
        end
        
    end
    
    cellMesh=[cellHeader,cellNodeCoord,cellConnCoord];
    
end

function [cellMesh]=CellCellCentredMesh(coordDat,cellCentredData,nodeConnDat,elmConnDat,vararg)
    % Saves the mesh as a cell-centred polygonal mesh
    
    [~,numNodes]=size(coordDat);
    [~,numElm]=size(cellCentredData);
    [~,numFaces]=size(elmConnDat);
    
    % Zone Header
    cellHeader=CellCentredDataHeader(numNodes,numElm,numFaces,vararg);
    % Coordinates
    cellStringNodal={};
    cellStringNodeConn={};
    cellStringElmConn={};
    if isempty(vararg.varsharezone)
        [cellStringNodal]=WriteNumArrayToCell(coordDat);
    end
    [cellStringCellCentred]=WriteNumArrayToCell(cellCentredData);
    if isempty(vararg.connecShareZone)
        [cellStringNodeConn]=WriteIntArrayToCell(nodeConnDat);
        [cellStringElmConn]=WriteIntArrayToCell(elmConnDat);
    end
    
    cellMesh=[cellHeader,cellStringNodal,cellStringCellCentred...
        ,cellStringNodeConn,cellStringElmConn];
    
end

function [cellStringArray]=WriteNumArrayToCell(array)
    
    sizArray=size(array);
    for ii=sizArray(1):-1:1
        cellStringArray{ii}='';
        for jj=1:sizArray(2)
            cellStringArray{ii}=...
                [cellStringArray{ii},num2str(array(ii,jj),24),'  '];
        end
        
    end
    
end

function [cellStringArray]=WriteIntArrayToCell(array)
    
    sizArray=size(array);
    for ii=sizArray(1):-1:1
        cellStringArray{ii}='';
        for jj=1:sizArray(2)
            cellStringArray{ii}=...
                [cellStringArray{ii},int2str(array(ii,jj)),'  '];
        end
        
    end
    
end

function []=ExtractCellCentredInfo(baseCellCentredGrid)
    
    for ii=1:length(baseCellCentredGrid)
        baseCellCentredGrid(ii).edge=[baseCellCentredGrid(ii).edge(:).index];
        baseCellCentredGrid(ii).numvertex=length(baseCellCentredGrid(ii).vertexindex);
    end
    
end


%% Various

function [cellCentredGrid]=CellCentredGridInformationReduced(refinedGrid)
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
    
    
end

