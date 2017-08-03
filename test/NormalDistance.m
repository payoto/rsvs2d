function [errorMeasure,areaDistrib]=NormalDistance(profileCoord,targCoord)
    
    
    normVec=@(v) v./repmat(sqrt(sum(v.^2,2)),[1,size(v,2)]);
    pDistDir=@(m,p,v) sign((m-repmat(p,[size(m,1),1]))*v').*sqrt(sum((m-repmat(p,[size(m,1),1])).^2,2));
    epsScales=0.01;
    
 
    maxDist=[min([profileCoord;targCoord]);max([profileCoord;targCoord])];
    nP1=size(profileCoord,1);
    edgeNormal=([0 1;-1 0]*normVec(profileCoord([2:end,1],:)-profileCoord)')';
    [~,vDist]=knnsearch(targCoord,profileCoord);
    for ii=1:nP1
        
        vertNormal=sum(edgeNormal(mod((ii-1:ii)-1,nP1)+1,:));
        
        scales=(maxDist-profileCoord([ii,ii],:))./repmat(vertNormal,[2,1]);
        scales=scales(:);
        scales=[min(scales(scales>0)),max(scales(scales<0))];
        
        if numel(scales)<2
            if isempty(scales)
                scales=[epsScales -epsScales];
            elseif scales<0
                scales=[epsScales,scales];
            elseif scales>0
                scales=[scales,-epsScales];
            end
        end
        testCoord=repmat(profileCoord(ii,:),[2,1])+(scales'*vertNormal);
        
        [x0,y0,iout,jout] = intersections(testCoord(:,1),testCoord(:,2),...
                targCoord([1:end,1],1),targCoord([1:end,1],2));
        rmv=isnan(iout) | isnan(jout);
        
        x0(rmv)=[];
        y0(rmv)=[];
        
        
        if numel(x0)>0
            vertDists=pDistDir([x0,y0],profileCoord(ii,:),vertNormal);
            
            vDist(ii)=min([vDist(ii);abs(vertDists)]);
        
        end
    end
    
    
    errorMeasure.sum=sum(vDist);
    errorMeasure.mean=mean(vDist);
    errorMeasure.std=std(vDist);
    errorMeasure.max=max(vDist);
    errorMeasure.min=min(vDist);
    areaDistrib=[profileCoord(:,1),vDist];
   
    
    %error('Compare Profile through area integration has not been coded yet')
end