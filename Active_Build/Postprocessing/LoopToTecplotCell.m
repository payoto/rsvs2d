function [cellWrite]=LoopToTecplotCell(loop,fieldpath)
   
    cellWrite=cell(0);
    for jj=1:numel(fieldpath)
        for ii=1:numel(loop)
            [loopdat]=ExtractFieldDat(loop(ii),fieldpath{jj});
            [cellMesh]=CoordToCellOut(loopdat);
            cellWrite=[cellWrite,cellMesh];
        end
    end
    
    
    
end


function [loop]=ExtractFieldDat(loop,fieldpath)
    
    for ii=1:numel(fieldpath)
        loop=loop.(fieldpath{ii});
    end
    
end

function [cellMesh]=CoordToCellOut(coord)
    
    
    coordDat=coord;
    vectorDat=coord;
    velDat=zeros([size(coord,1),1]);
    %velDat=[[snaxel(:).isfreeze]']*2+1;
    %velDat=[[snaxel(:).orderedge]'];
    vectorDat(:,1)=vectorDat(:,1).*velDat;
    vectorDat(:,2)=vectorDat(:,2).*velDat;
    vectorDat=[vectorDat,velDat];
    vertIndex=1:size(coord,1);
    connDat=[[1:size(coord,1)]',(mod(1:size(coord,1),size(coord,1))+1)'];
    
    [cellMesh]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,1);
    
end