function [] = include_GridCheck()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


function []=CheckGrid(gridreshape)
    domainBounds=[0 1 ; 0 1];
    
    [unstructured]=ModifReshape(gridreshape);
    
    
    figh=figure;
    axh=axes;
    hold on
    
    colString='bgcmyk';
    
    isEdgeSub=find(unstructured.edge.index);
    for ii=1:length(isEdgeSub)
        PlotEdgeGrid(figh,axh,unstructured,isEdgeSub(ii),'b-')
    end
    
%     isEdgeSub=find(~unstructured.edge.boundaryis0);
%     for ii=1:length(isEdgeSub)
%         PlotEdgeGrid(figh,axh,unstructured,isEdgeSub(ii),'b-')
%     end
    
    
    axis equal
    axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    
    
end

function []=PlotEdgeGrid(figh,axh,unstructured,subEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(subEdge,:);
    vertsub(1)=find(unstructured.vertex.index==vertices(1));
    vertsub(2)=find(unstructured.vertex.index==vertices(2));
    coord=unstructured.vertex.coord(vertsub,:);
    
    plot(coord(:,1),coord(:,2),format)
    text(mean(coord(:,1)),mean(coord(:,2)),int2str(unstructured.edge.index(subEdge)),'color','b')
    text(mean(coord(1,1)),mean(coord(1,2)),int2str(vertices(1)),'color','g')
    text(mean(coord(2,1)),mean(coord(2,2)),int2str(vertices(2)),'color','g')
end

function []=PlotLoopGrid(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).vertex.coord;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotSubDivGrid(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).subdivision;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotCellGrid(figh,axh,unstructured,indexCell,format)
    figure(figh)
    axes(axh)
    
    
    coord=unstructured.cell.coord(indexCell,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end

