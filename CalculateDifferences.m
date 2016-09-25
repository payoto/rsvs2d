

function []=CalculateDifferences(sensloops,snakloops)
    
    if numel(sensloops)~=numel(snakloops)
        error('Size mismatch beween the two sources')
    end
    
    figure,
    
    
    for ii=1:numel(sensloops)
        
        errstruct(ii).coordsens=sensloops(ii).loop.subdivspline;
        errstruct(ii).coordsnak=snakloops(ii).loop.subdivspline;
        errstruct(ii).diff=errstruct(ii).coordsens-errstruct(ii).coordsnak;
        errstruct(ii).dist=sqrt(sum(errstruct(ii).diff.^2,2));
        errstruct(ii).sum=sum(errstruct(ii).dist);
        errstruct(ii).mean=mean(errstruct(ii).dist);
        errstruct(ii).std=std(errstruct(ii).dist);
        errstruct(ii).min=min(errstruct(ii).dist);
        errstruct(ii).max=max(errstruct(ii).dist);
        
        semilogy(errstruct(ii).coordsnak(:,1),errstruct(ii).dist)
        hold on
    end
    figure
    for ii=1:numel(sensloops)
        semilogy(errstruct(ii).coordsnak(:,2),errstruct(ii).dist)
        hold on
    end
    figure
    for ii=1:numel(sensloops)
        semilogy(1:numel(errstruct(ii).dist),errstruct(ii).dist)
        hold on
    end
    figure,
    semilogy([errstruct(:).sum])
    figure,
    semilogy([errstruct(:).mean])
    ylabel('mean distance')
    figure,
    semilogy([errstruct(:).std])
    ylabel('std distance')
    figure,
    semilogy([errstruct(:).min])
    ylabel('min distance')
    figure,
    semilogy([errstruct(:).max])
    ylabel('max distance')
    
    
    
end