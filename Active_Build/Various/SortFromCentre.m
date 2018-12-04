function [ind]=SortFromCentre(numin)
    flagReshape=0;
    if(size(numin,1)>1)
        numin=numin';
        flagReshape=1;
    end
    if mod(numel(numin),2)==0
        [ind]=SortFromCentreEven(numin);
    else
        [ind]=SortFromCentreOdd(numin);
    end
    
    if flagReshape
        ind=ind';
    end
end


function [ind]=SortFromCentreOdd(numin)
    
    [~,i]=max(numin);
    ex=0;
    if i>numel(numin)/2
        ex=1;
    end
    midpoint = ceil(numel(numin)/2)-ex;
    [~,indStart]=sort(numin(1:midpoint));
    [~,indEnd]=sort(numin(midpoint+1:end));
    
    ind = [indStart,flip(indEnd)+midpoint];
    
end

function [ind]=SortFromCentreEven(numin)
    

    
    [~,indStart]=sort(numin(1:numel(numin)/2));
    [~,indEnd]=sort(numin(numel(numin)/2+1:end));
    
    ind = [indStart,flip(indEnd)+numel(numin)/2];
    
end