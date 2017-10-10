% THis code has been integrated in GridRefinement.m

function [refGrid]=RefineTriangular(baseGrid,nRefine)
    
    tempGrid=baseGrid;
    for ii=1:nRefine
        
       [tempGrid,cellVert]=BreakEdgesHalf(tempGrid);
       [tempGrid]=AddNewEdgesCell(tempGrid,cellVert);
    end
    
    
    refGrid=tempGrid;
end




function [baseGrid,cellVert]=BreakEdgesHalf(baseGrid)
    
    vertexInd=[baseGrid.vertex.index];
    vertexCoord=vertcat(baseGrid.vertex.coord);
    cellInd=[0,baseGrid.cell.index];
    
    maxVertInd=max(vertexInd)+1;
    maxEdgeInd=max([baseGrid.edge.index])+1;
    cellVert.cell=repmat(struct('index',[],'vertex',[]),[numel(baseGrid.cell)+1,1]);
    cellVert.vertex=repmat(struct('index',[],'edge',[]),[numel(baseGrid.edge),1]);
    
    for ii=1:numel(baseGrid.edge)
        
        vertSub=FindObjNum([],baseGrid.edge(ii).vertexindex,vertexInd);
        cellSub=FindObjNum([],baseGrid.edge(ii).cellindex,cellInd);
        
        kk=numel(baseGrid.edge)+1;
        baseGrid.edge(kk)=baseGrid.edge(ii);
        baseGrid.edge(kk).index=maxEdgeInd;
        
        kkv=numel(baseGrid.vertex)+1;
        baseGrid.vertex(kkv).index=maxVertInd;
        
        baseGrid.vertex(kkv).coord=sum(vertexCoord(vertSub,:))/2;
        
         baseGrid.edge(kk).vertexindex(1)=maxVertInd;
         baseGrid.edge(ii).vertexindex(2)=maxVertInd;
         for jj=1:numel(cellSub)
             cellVert.cell(cellSub(jj)).index=cellInd(cellSub(jj));
             cellVert.cell(cellSub(jj)).vertex=[cellVert.cell(cellSub(jj)).vertex,maxVertInd];
         end
         cellVert.vertex(ii).index=maxVertInd;
         cellVert.vertex(ii).edge=[baseGrid.edge(ii).index,maxEdgeInd];
         
         maxVertInd=maxVertInd+1;
         maxEdgeInd=maxEdgeInd+1;
    end

end

function [refGrid]=AddNewEdgesCell(baseGrid,cellVert)
    
    cellInd=[0,baseGrid.cell.index];
    edgeInd=[baseGrid.edge.index];
    
    cellEdgeInd=[baseGrid.edge.cellindex];
    vertEdgeInd=[baseGrid.edge.vertexindex];
    
    maxEdgeInd=max(edgeInd)+1;
    maxCellInd=max(cellInd)+1;
    
    vertCellVert=[cellVert.vertex.index];
    vertCellCell=[cellVert.cell.index];
    nEdge=numel(baseGrid.edge);
    
    baseGrid.edge(end+3*numel(baseGrid.cell)).index=[];
    for ii=1:numel(baseGrid.cell)
        
        subCellVert=FindObjNum([],baseGrid.cell(ii).index,vertCellCell);
        newCellInd=[maxCellInd:maxCellInd+2];
        
        newEdgeInd=[maxEdgeInd:maxEdgeInd+2]';
        edgeCellIndNew=[ones(3,1)*baseGrid.cell(ii).index,newCellInd'];
        edgeVertIndNew=cellVert.cell(subCellVert).vertex([1 2;2 3;3 1]);
        
        newEdgeSub=nEdge:nEdge+2;
        for jj=1:3;
            nEdge=nEdge+1;
            baseGrid.edge(nEdge).index=newEdgeInd(jj);
            baseGrid.edge(nEdge).vertexindex=edgeVertIndNew(jj,:);
            baseGrid.edge(nEdge).cellindex=edgeCellIndNew(jj,:);
        end
        
        % for each new Edge need to find the two old edge which closes the
        % containment and create the new cell
        % these are the ones which share a vertex
        
        %edgeCellSub=ceil(FindObjNum([],baseGrid.cell(ii).index,cellEdgeInd)/2);
        for jj=1:3
            vertSub=FindObjNum([],edgeVertIndNew(jj,:),vertCellVert);
            edgeSub=FindObjNum([],[cellVert.vertex(vertSub).edge],edgeInd);
            vertEdgeSub=[baseGrid.edge(edgeSub).vertexindex];
            vertEdgeSub=vertEdgeSub((vertEdgeSub~=edgeVertIndNew(jj,1)) & ...
                (vertEdgeSub~=edgeVertIndNew(jj,2)));
            cellSim=FindIdenticalVector(vertEdgeSub');
            edgeSub=edgeSub(cellSim{cellfun(@numel,cellSim)>1});
            for kk=1:2
                baseGrid.edge(edgeSub(kk)).cellindex(...
                    baseGrid.edge(edgeSub(kk)).cellindex==...
                    baseGrid.cell(ii).index)=newCellInd(jj);
            end
            baseGrid.cell(end+1)=baseGrid.cell(ii);
            baseGrid.cell(end).index=newCellInd(jj);
        end
        
        % And then create the new cells
        
        maxCellInd=maxCellInd+3;
        maxEdgeInd=maxEdgeInd+3;
    end
    refGrid=baseGrid;
end
