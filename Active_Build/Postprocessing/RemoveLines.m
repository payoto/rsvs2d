function []=RemoveLines(figh,lStyle)
    
    
    
    for ii=1:numel(figh)
        
        l=findobj(figh(ii),'type','line');
        isDel=false(size(l));
        for jj=1:numel(l)
            isDel(jj)=strcmp(l(jj).LineStyle,lStyle);
        end
        delete(l(isDel))
    end
    
    
    
end