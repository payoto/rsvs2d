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

%#codegen

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
        volumefraction(ii).isSnax=false;
        jj=1;
        nSearch=length(newSubs);
        while (~volumefraction(ii).isSnax && jj<=nSearch)
            
            volumefraction(ii).isSnax=volumefraction(ii).isSnax  ...
                || ~isempty(cellCentredGrid(newSubs(jj)).snaxel);
            
            jj=jj+1;
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
        snaxCells(snaxCells==0)=[];
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
    
    actCoord=zeros(size([2,length(snakSub)]));
    
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
