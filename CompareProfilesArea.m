function [errorMeasure,areaDistrib]=CompareProfilesArea(profileCoord,targCoord)
    
    [x0,y0,iout,jout] = intersections(profileCoord([1:end,1],1),profileCoord([1:end,1],2),...
        targCoord([1:end,1],1),targCoord([1:end,1],2));
    rmv=isnan(iout) | isnan(jout);
    iout(rmv)=[];
    x0(rmv)=[];
    y0(rmv)=[];
    jout(rmv)=[];
    
    [iout,sortIndex]=sort(iout);
    x0=x0(sortIndex);
    y0=y0(sortIndex);
    jout=jout(sortIndex);
    targArea=abs(CalculatePolyArea(targCoord));
    nX0=numel(x0);
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    
    if nX0>0
        areaErr=zeros(size(x0))';
        areaLength=zeros(size(x0))';
        areaPosx=zeros(size(x0))';
        areaPosXmin=zeros(size(x0))';
        areaPosXmax=zeros(size(x0))';
        for ii=1:nX0
            % Need to make sure both profiles go in the same direction
            iip1=mod(ii,nX0)+1;
            if floor(iout(ii))~=floor(iout(iip1))
                if ceil(iout(ii))>floor(iout(iip1))
                    ind1=[ceil(iout(ii)):nP1,1:floor(iout(iip1))];
                else
                    ind1=ceil(iout(ii)):floor(iout(iip1));
                end
            else
                ind1=[];
            end
            if floor(jout(ii))~=floor(jout(iip1))
                if ceil(jout(ii))>floor(jout(iip1))
                    ind2=flip([ceil(jout(ii)):nP2,1:floor(jout(iip1))]);
                else
                    ind2=flip(ceil(jout(ii)):floor(jout(iip1)));
                end
            else
                ind2=[];
            end
            
            actPts=[[x0(ii),y0(ii)];
                profileCoord(ind1,:);
                [x0(iip1),y0(iip1)];
                targCoord(ind2,:)];
            areaErr(ii)=abs(CalculatePolyArea(actPts));
            areaLength(ii)=max(actPts(:,1))-min(actPts(:,1));
            areaPosx(ii)=(min(actPts(:,1))+max(actPts(:,1)))/2;
            areaPosy(ii)=(min(actPts(:,2))+max(actPts(:,2)))/2;
            areaPosXmin(ii)=min(actPts(:,1));
            areaPosXmax(ii)=max(actPts(:,1));
        end
        
        errorMeasure.sum=sum(areaErr)/targArea;
        errorMeasure.mean=mean(areaErr)/targArea;
        errorMeasure.std=std(areaErr)/targArea;
        errorMeasure.max=max(areaErr)/targArea;
        errorMeasure.min=min(areaErr)/targArea;
        
        areaDistrib=[min(areaPosXmin),areaPosx,areaPosXmax;
            0,areaErr./areaLength*2,-areaErr./areaLength*2];
        [areaDistrib(1,:),iSortDistrib]=sort(areaDistrib(1,:));
        eqInd=find(abs(areaDistrib(1,1:end-1)-areaDistrib(1,2:end))<1e-10);
        areaDistrib(2,:)=(areaDistrib(2,iSortDistrib));
        areaDistrib(2,eqInd)=areaDistrib(2,eqInd)+areaDistrib(2,eqInd+1);
        areaDistrib(:,eqInd+1)=[];
        areaDistrib(2,:)=cumsum(areaDistrib(2,:));
        areaDistrib(:,find(areaDistrib(2,:)<1e-18))=[];
        areaDistrib=areaDistrib';
        %         distribError=areaErr./areaLength;
        %         xPts=[[profileCoord(1,1),reshape(x0,[1,numel(x0)])];
        %             [reshape(x0,[1,numel(x0)]),profileCoord(end,1)]];
        %         xPts=xPts(:);
        %         errPts=[distribError(end),distribError;distribError(end),distribError];
        %         errPts=errPts(:);
        %         errorDist=[xPts',errPts'];
    else
        
        [errorMeasure]=(abs(CalculatePolyArea(profileCoord))-targArea)/targArea;
        areaDistrib=[(min(profileCoord(:,1))+min(profileCoord(:,1)))/2,errorMeasure];
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end