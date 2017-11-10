function [refinedGrid]=MakeTriangles(refinedGrid)
    % function adds diagonal elements to a rectangular grid.
    
    vertInd=[refinedGrid.vertex(:).index];
    edgeInd=[refinedGrid.edge(:).index];
    cellInd=[refinedGrid.cell(:).index];
    cellCentredGrid=CellCentredGrid(refinedGrid);
    
    maxEdge=max(edgeInd)+1;
    maxCell=max(cellInd)+1;
    for ii=1:numel(cellCentredGrid)
        cellVert=[cellCentredGrid.vertex(:).index];
        cellEdge=[cellCentredGrid.edge(:).index];
        edgeVert=[cellCentredGrid.edge(:).vertexindex];
        vert1=cellVert(1);
        
        edgeVert(1:2:end)==vert1 | edgeVert(2:2:end)==vert1
        
        
    end
    
    
end

