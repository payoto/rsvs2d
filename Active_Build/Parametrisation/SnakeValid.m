function [out]=SnakeValid()
    
    
    snakCell={'Snakestestsmooth1','Snakestestsmooth1_2','Snakestestsmooth2',...
        'Snakestestsmooth3','Snakestestsmooth3_1','Donught','Donught2',...
        'SnakesFoilVVSmall','BuzmanBiplane3','SnakesFoilVVSmall4',...
        'WeirdShapeIn','WeirdShapeOut'};
    
    T{length(snakCell)}=[];
    
    out=repmat(struct('unstructured',[],'loop',[],'unstructReshape',[],'snakSave',[]...
    ),[1 length(snakCell)]);
    
    for ii=1:length(snakCell)
        [T{ii},out(ii)]=CallMain(snakCell{ii});
        
        T{ii}
    end
    
    
    
    
    
    
end

function [T,out]=CallMain(snakCell)
    
    [T,out.unstructured,out.loop,out.unstructReshape,out.snakSave]=...
    evalc('Main([''val_'',snakCell])');
end