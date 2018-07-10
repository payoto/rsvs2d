function []=CompareASOLoops(ASOstruct,iASostruct)
    
    
    plotPoints= @(points,c) plot(points([1:end],1),points([1:end],2),'color',c);
    plotLines=@(vec1,vec2) plot([vec1(:,1),vec2(:,1)]',[vec1(:,2),vec2(:,2)]','k--');
    
    for kk=1:numel(iASostruct)
        figure, hold on
        c=get(gca,'colororder');
        maxPerLoop=zeros(size(ASOstruct(iASostruct(kk)).loops{1}));
        for ii=1:numel(ASOstruct(iASostruct(kk)).loops{1})
            [~,min1]=min(ASOstruct(iASostruct(kk)).loops{1}{ii}(:,1));
            [~,min2]=min(ASOstruct(iASostruct(kk)).loops{end}{ii}(:,1));
            loopStart=ASOstruct(iASostruct(kk)).loops{1}{ii}([min1:end,1:min1-1],:);
            loopEnd=ASOstruct(iASostruct(kk)).loops{end}{ii}([min2:end,1:min2-1],:);
            
            plotLines(loopEnd,loopStart)
            plotPoints(loopStart,c(1,:))
            plotPoints(loopEnd,c(2,:))
            vec=loopStart-loopEnd;
            normVec=sqrt(sum(vec.^2,2));
            
            maxPerLoop(ii)=max(normVec);
        end
        maxPerLoop
    end
end