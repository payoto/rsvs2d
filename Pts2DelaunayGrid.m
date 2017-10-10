% This code has been adapted into GridGeneration

function [gridtri]=Pts2DelaunayGrid(pts)
    
    dtri=delaunayTriangulation(pts);
    [gridtri]=BasicGrid(dtri,pts);
    
    %CheckGrid(gridtri)
end



function [gridtri]=BasicGrid(dtri,pts)
    
    kke=1;
    ndtri2=size(dtri,2);
    gridtri.cell=repmat(struct('index',0,'fill',[]),[size(dtri,1),1]);
    gridtri.vertex=repmat(struct('index',0,'coord',[0 0]),[size(pts,1),1]);
    gridtri.edge=repmat(struct('index',0,'cellindex',[],'vertexindex',[]),[size(dtri,1)*3,1]);
    for ii=1:size(pts,1)
        gridtri.vertex(ii).index=ii;
        gridtri.vertex(ii).coord=pts(ii,:);
    end
    for ii=1:size(dtri,1)
        gridtri.cell(ii).index=ii;
        for jj=1:ndtri2
            gridtri.edge(kke).vertexindex=sort([dtri(ii,jj),dtri(ii,mod(jj,ndtri2)+1)]);
            gridtri.edge(kke).cellindex=ii;
            
            kke=kke+1;
        end
    end
    
    
    cellSameEdge=FindIdenticalVector(vertcat(gridtri.edge(:).vertexindex));
    edgeRm=[];
    for ii=1:numel(cellSameEdge)
        gridtri.edge(cellSameEdge{ii}(1)).cellindex=...
            [gridtri.edge(cellSameEdge{ii}).cellindex];
        if numel(gridtri.edge(cellSameEdge{ii}(1)).cellindex)==1
           gridtri.edge(cellSameEdge{ii}(1)).cellindex(2)=0; 
        end
        gridtri.edge(cellSameEdge{ii}(1)).index=ii;
        edgeRm=[edgeRm,cellSameEdge{ii}(2:end)];
    end
    gridtri.edge(edgeRm)=[];
    
end