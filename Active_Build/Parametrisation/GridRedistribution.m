%% Grid redistributions
function [unstructReshape,gridrefined,connectstructinfo,loop,...
        unstructured, unstructuredrefined]=GridRedistribution(...
        gridrefined,connectstructinfo,unstructReshape,loop,param)
        
    varExtract={'gridDistrib', 'randMultiplier', 'boundstr'};
    [gridDistrib, randMultiplier, boundstr]=ExtractVariables(varExtract,param);
    exec = true;
    switch gridDistrib
        case 'randMod'
            [unstructReshape,gridrefined,connectstructinfo]=TestNonRectangleGrid...
                (gridrefined,connectstructinfo,unstructReshape, randMultiplier);
        case 'waveMod'
            [unstructReshape,gridrefined,connectstructinfo]=WaveNonRectangleGrid(gridrefined,...
                connectstructinfo, unstructReshape);
        otherwise
            exec = false;
    end
    
    
    if exec
        [gridrefined]=EdgePropertiesReshape(gridrefined);
        [loop]=GenerateSnakStartLoop(gridrefined,boundstr);
    end
    [unstructured]=ModifReshape(unstructReshape);  
    [unstructuredrefined]=ModifReshape(gridrefined);
end
function [grid] = ApplyCoordTransform(grid, coordtrans)
    
    for ii = 1:numel(grid.vertex)
        grid.vertex(ii).coord =...
            coordtrans(grid.vertex(ii).coord);

    end
end

function [gridbase,gridrefined,connectstructinfo]=WaveNonRectangleGrid(gridrefined,...
    connectstructinfo, gridbase)
    
    [gridbase,coarsenconnec]=CoarsenGrid(gridrefined,gridbase,connectstructinfo);
    gridbase.vertex=gridbase.vertex(FindObjNum([],...
        unique([gridbase.edge.vertexindex]),[gridbase.vertex.index]));
    coordtrans= @(coord) coord + [zeros(size(coord,1),1),...
                sin(coord(:,1)*3*pi)*0.1];
    [gridrefined] = ApplyCoordTransform(gridrefined, coordtrans);
    [gridbase.vertex.coord]=deal(gridrefined.vertex(...
        FindObjNum([],[gridbase.vertex.index],[gridrefined.vertex.index])...
        ).coord);
    
end

function [gridbase,gridrefined,connectstructinfo]=TestNonRectangleGrid(gridrefined,...
    connectstructinfo, gridbase, randMultiplier)
    
    [gridbase,coarsenconnec]=CoarsenGrid(gridrefined,gridbase,connectstructinfo);
    gridbase.vertex=gridbase.vertex(FindObjNum([],...
        unique([gridbase.edge.vertexindex]),[gridbase.vertex.index]));
    
    [gridrefined]=RandGridDistrib(gridrefined, randMultiplier);
    [gridbase.vertex.coord]=deal(gridrefined.vertex(...
        FindObjNum([],[gridbase.vertex.index],[gridrefined.vertex.index])...
        ).coord);
    
end


function [unstructured]=RandGridDistrib(unstructured, randMultiplier)
    
    coord=vertcat(unstructured.vertex.coord);
    
    for ii=1:size(coord,2)
        x=coord(:,ii);
        x=unique(x);
        Dx=x(2:end)-x(1:end-1);
        Dx=min(Dx(Dx>1e-10));
        coord(:,ii)=Dx*rand(size(coord(:,1)))*randMultiplier/2+coord(:,ii);
    end
    for ii=1:size(coord,1)
        unstructured.vertex(ii).coord=coord(ii,:);
    end
    
end