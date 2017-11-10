

function []=CalculateDifferences(sensloops,snakloops)
    % Calculates differences between 2 loops
    if numel(sensloops)~=numel(snakloops)
        error('Size mismatch beween the two sources')
    end
    plotPoints= @(points,str) plot(points([1:end],1),points([1:end],2));
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
    end
    
%     figure,
%     for ii=1:numel(sensloops)   
%         semilogy(errstruct(ii).coordsnak(:,1),errstruct(ii).dist)
%         axis([0 1 1e-12 1])
%         title(sprintf('Mode %i',(ii-2-mod(ii-2,2))/2+1))
%         pause; 
%     end
%     ylabel('distance')
%     xlabel('x')
    figure,
    for ii=1:numel(sensloops)   
        semilogy(errstruct(ii).coordsnak(:,1),errstruct(ii).dist)
        hold on
    end
    ylabel('distance')
    xlabel('x')
    figure
    for ii=1:numel(sensloops)
        semilogy(errstruct(ii).coordsnak(:,2),errstruct(ii).dist)
        hold on
    end
    ylabel('distance')
    xlabel('y')
    figure
    for ii=1:numel(sensloops)
        semilogy(1:numel(errstruct(ii).dist),errstruct(ii).dist)
        hold on
    end
    ylabel('distance')
    xlabel('index')
    figure,
    semilogy([errstruct(:).sum])
    ylabel('sum distance')
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
    eps=1e-6;
    plotProfiles=find([errstruct(:).max]>eps);
    for ii=plotProfiles
        figure,
        hold on
        plotPoints(errstruct(1).coordsnak)
        plotPoints(errstruct(ii).coordsens)
        plotPoints(errstruct(ii).coordsnak)
        [maxInd]=find(errstruct(ii).dist>eps);
        plot(errstruct(ii).coordsens(maxInd,1),errstruct(ii).coordsens(maxInd,2),'r*')
        plot(errstruct(ii).coordsnak(maxInd,1),errstruct(ii).coordsnak(maxInd,2),'k*')
        title(num2str(errstruct(ii).max))
        pause
        close
    end
        
    
end