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


function [cellCentredGrid]=CellCentreGridInformation(refinedGrid)
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
        orderedVertexSub=FindObjNum([],orderedVertices,[cellCentredGrid(ii).vertex(:).index]);
        coord=vertcat(cellCentredGrid(ii).vertex(orderedVertexSub).coord);
        [isCCW]=CCWLoop(coord);
        for jj=1:length(orderedEdges)
            
            edgeVertSub=FindObjNum([],cellOrderedVertex{1}(jj,:),[cellCentredGrid(ii).vertex(:).index]);
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
