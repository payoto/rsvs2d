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
    cellVert.cell=repmat(struct(vertex,[]),[numel(baseGrid.cell)+1,1]);
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
    maxEdgeInd=max([baseGrid.edge.index])+1;
    maxCellInd=max(cellInd)+1;
    
    
end
