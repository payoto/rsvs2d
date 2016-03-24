
function [] = include_EdgeInformation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%% Reshaped edge information

function [loop]=OrderSurfaceVertexReshape(gridreshape,isEdge,cond)
    % function ordering the surface vertices in counter-clockwise order
    % isEdge is a logical indexing arraying informing which edges are on the
    % edge of a surface
    % cond is the optional string argument describing the condition fulfilled
    % by isEdge possible arg: '0bound', '1bound', 'intermBound'
    
    isEdgeIndex=find(isEdge);
    
    blockEdges=vertcat(gridreshape.edge(isEdgeIndex).vertexindex);
    blockCell=vertcat(gridreshape.edge(isEdgeIndex).cellindex);
    fillCell=vertcat(gridreshape.edge(isEdgeIndex).fill);
    if ~exist('cond','var'), cond='1bound';end
    switch cond
        case '0bound'
            fillCell=fillCell>0;
        case '1bound'
            fillCell=fillCell==1;
        case 'intermBound'
            fillCell(:,1)=fillCell(:,1)>fillCell(:,2);
            fillCell(:,2)=fillCell(:,2)>fillCell(:,1); 
    end
    for ii=1:length(fillCell(:,1))
        %colNum=find(fillCell(ii,:));
        blockCellTrunc(ii)=blockCell(ii,find(fillCell(ii,:)));
        
    end
    
    coordVertex=vertcat(gridreshape.vertex(:).coord);
    vertexIndex=[gridreshape.vertex(:).index];
    % Order edges into closed loops
    %[cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(blockEdges,blockCellTrunc);
    [cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(blockEdges);
    
    for ii=1:length(cellOrderedVertex)
        loop(ii).vertex.index=[cellOrderedVertex{ii}(:,1);cellOrderedVertex{ii}(1:2,1)];
        loopVertSub=FindObjNum([],loop(ii).vertex.index,vertexIndex);
        loop(ii).vertex.coord=coordVertex(loopVertSub,:);
        loop(ii).edge.index=[gridreshape.edge(isEdgeIndex(cellOrderedEdges{ii})).index];
    end
    
end

function [gridreshape]=EdgePropertiesReshape(gridreshape)
    % Extracts fill data for each edge and classifies edges depending on
    % neighbouring cells
    
    % Extract fill information
    edgeFill=EdgeFillInformationReshape(gridreshape);
    
    % Extract indices of edges matching boundary criteria
    edgeFillSort=sort(edgeFill,2);
    
    is0=sum((edgeFillSort==0),2);
    is1=sum((edgeFillSort==1),2);
    isElse=sum((edgeFillSort~=0)&(edgeFillSort~=1),2);
    
    boundaryisHalf=((isElse==2));
    boundaryis0=((is0==1));
    boundaryis1=((is1==1));

    solidisIn0=((is0==2));
    solidnotIn0=((is0~=2));
    solidisIn1=((is1==2));
    
    for ii=1:length(gridreshape.edge)
        [gridreshape.edge(ii).boundaryisHalf]=boundaryisHalf(ii);
        [gridreshape.edge(ii).boundaryis0]=boundaryis0(ii);
        [gridreshape.edge(ii).boundaryis1]=boundaryis1(ii);

        [gridreshape.edge(ii).solidisIn0]=solidisIn0(ii);
        [gridreshape.edge(ii).solidnotIn0]=solidnotIn0(ii);
        [gridreshape.edge(ii).solidisIn1]=solidisIn1(ii);
    
        [gridreshape.edge(ii).fill]=edgeFill(ii,:);
    end
%     unstructured.edge.boundary=boundary;
%     unstructured.edge.solid=solid;
    
end

function [edgeFill]=EdgeFillInformationReshape(gridreshape)
    % returns the fill information of neighbouring cells
    
    
    % decomposing structure
    fillCellDat=[0;vertcat(gridreshape.cell(:).fill)];
    edgeCellIndex=vertcat(gridreshape.edge(:).cellindex);
    cellIndex=[gridreshape.cell(:).index];
    
    edgeCellSub=FindObjNum(gridreshape.cell,edgeCellIndex(:,1),cellIndex);
    edgeCellSub(:,2)=FindObjNum(gridreshape.cell,edgeCellIndex(:,2),cellIndex);
    % detecting array size
    [mCI,nCI]=size(edgeCellIndex);
    
    % Preallocating
    edgeFill=zeros(mCI,nCI);
    
    for ii=1:nCI
        edgeFill(:,ii)=fillCellDat(edgeCellSub(:,ii)+1);
    end
    
end

%% Surface Identification functions

function [loop]=OrderSurfaceVertex(unstructured,isEdge,cond)
    % function ordering the surface vertices in counter-clockwise order
    % isEdge is a logical indexing arraying informing which edges are on the
    % edge of a surface
    % cond is the optional string argument describing the condition fulfilled
    % by isEdge possible arg: '0bound', '1bound', 'intermBound'
    
    isEdgeIndex=find(isEdge);
    
    blockEdges=unstructured.edge.vertexindex(isEdgeIndex,:);
    blockCell=unstructured.edge.cellindex(isEdgeIndex,:);
    fillCell=unstructured.edge.fill(isEdgeIndex,:);
    if ~exist('cond','var'), cond='1bound';end
    switch cond
        case '0bound'
            fillCell=fillCell>0;
        case '1bound'
            fillCell=fillCell==1;
        case 'intermBound'
            fillCell(:,1)=fillCell(:,1)>fillCell(:,2);
            fillCell(:,2)=fillCell(:,2)>fillCell(:,1); 
    end
    for ii=1:length(fillCell(:,1))
        %colNum=find(fillCell(ii,:));
        blockCellTrunc(ii)=blockCell(ii,find(fillCell(ii,:)));
        
    end
    
    coordVertex=unstructured.vertex.coord;
    
    % Order edges into closed loops
    %[cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(blockEdges,blockCellTrunc);
    [cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(blockEdges);
    
    for ii=1:length(cellOrderedVertex)
        loop(ii).vertex.index=[cellOrderedVertex{ii}(:,1);cellOrderedVertex{ii}(1:2,1)];
        loop(ii).vertex.coord=coordVertex(loop(ii).vertex.index,:);
        loop(ii).edge.index=isEdgeIndex(cellOrderedEdges{ii});
    end
    
end

function [unstructured]=EdgeProperties(unstructured)
    % Extracts fill data for each edge and classifies edges depending on
    % neighbouring cells
    
    % Extract fill information
    edgeFill=EdgeFillInformation(unstructured);
    
    % Extract indices of edges matching boundary criteria
    edgeFillSort=sort(edgeFill,2);
    
    is0=sum((edgeFillSort==0),2);
    is1=sum((edgeFillSort==1),2);
    isElse=sum((edgeFillSort~=0)&(edgeFillSort~=1),2);
    
    unstructured.edge.boundaryisHalf=(isElse==2);
    unstructured.edge.boundaryis0=(is0==1);
    unstructured.edge.boundaryis1=(is1==1);
    
    unstructured.edge.solidisIn0=(is0==2);
    unstructured.edge.solidnotIn0=(is0~=2);
    unstructured.edge.solidisIn1=(is1==2);
    
    unstructured.edge.fill=edgeFill;
%     unstructured.edge.boundary=boundary;
%     unstructured.edge.solid=solid;
    
end

function [edgeFill]=EdgeFillInformation(unstructured)
    % returns the fill information of neighbouring cells
    
    
    % decomposing structure
    fillCellDat=[0;unstructured.cell.fill];
    edgeCellIndex=unstructured.edge.cellindex+1;
    
    % detecting array size
    [mCI,nCI]=size(edgeCellIndex);
    
    % Preallocating
    edgeFill=zeros(mCI,nCI);
    
    for ii=1:nCI
        edgeFill(:,ii)=fillCellDat(edgeCellIndex(:,ii));
    end
    
end

