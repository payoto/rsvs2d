function [fill,constrVal]=LoopToFill(loop,gridReshape)
    
    cellCentredGrid=CellCentreGridInformation(gridReshape);
    [cellCentredGrid]=CellPolygon(cellCentredGrid);
    
    fillAll=zeros([1,numel(cellCentredGrid)]);
    
    for ii=1:numel(loop)
        for jj=1:numel(cellCentredGrid)
            [internLoop]=ANDProfileOP(loop(ii).coord,cellCentredGrid(jj).ordered.coord);
            for kk=1:numel(internLoop)
                fillAll(jj)=fillAll(jj)+abs(CalculatePolyArea(internLoop(kk).coord));
            end
        end
    end
    fillAll=fillAll./[cellCentredGrid(:).volume];
    fill=fillAll(logical([cellCentredGrid(:).isactive]));
    
    constrVal={find(fill~=0),fill(fill~=0)};
end

function [cellGrid]=CellPolygon(cellGrid)
    
    for ii=1:numel(cellGrid)
        
        cellVertInd=vertcat(cellGrid(ii).edge(:).vertexindex);
        cellVertIndOrd=zeros(size(cellVertInd));
        cellEdgeOrd=zeros([size(cellVertInd,1),1]);
        kk=1;
        ll=1;
        for jj=1:size(cellVertInd,1)
            cellVertIndOrd(jj,ll:(3-2*ll):(3-ll))=cellVertInd(kk,:);
            cellEdgeOrd(jj)=cellGrid(ii).edge(kk).index;
            kkm1=kk;
            [kk,ll]=find(cellVertInd(kk,3-ll)==cellVertInd([1:kk-1,kk+1:end],:));
            kk=kk+(kk>=kkm1);
        end
        cellGrid(ii).ordered.edge=cellEdgeOrd;
        cellGrid(ii).ordered.vertex=cellVertIndOrd(:,1);
        cellGrid(ii).ordered.coord=...
            vertcat(cellGrid(ii).vertex(FindObjNum([],cellGrid(ii).ordered.vertex,...
            [cellGrid(ii).vertex(:).index])).coord);
    end
    
end

function [newloop]=ANDProfileOP(profileCoord,targCoord)
    % compares two loops and returns the area filled by both
    if ~CCWLoop(profileCoord)
        profileCoord=flip(profileCoord);
    end
    if ~CCWLoop(targCoord)
        targCoord=flip(targCoord);
    end
    
    [x0,y0,iout,jout] = intersections(profileCoord([1:end,1],1),profileCoord([1:end,1],2),...
         targCoord([1:end,1],1),targCoord([1:end,1],2));
%     [x0,y0,iout,jout] = intersections(profileCoord([1:end],1),profileCoord([1:end],2),...
%         targCoord([1:end],1),targCoord([1:end],2));
    rmv=isnan(iout) | isnan(jout);
    iout(rmv)=[];
    x0(rmv)=[];
    y0(rmv)=[];
    jout(rmv)=[];
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    iout=mod(iout-1,nP1)+1;
    jout=mod(jout-1,nP2)+1;
    
    [iout,sortIndex]=unique(iout);
    x0=x0(sortIndex);
    y0=y0(sortIndex);
    jout=jout(sortIndex);
    
    nX0=numel(x0);
    
    interCoord=[x0,y0];
    if nX0>0
        newloop=repmat(struct('coord',zeros(0,2)),[1,nX0]);
        nXplore=1:nX0;
        pLoop=0;
        while ~isempty(nXplore)
            % Need to make sure both profiles go in the same direction
            
            % THis code is insufficient it does not check for that no other
            % intersection lies on ind2 (cannot lie on ind1 as it is sorted)
            kk=1;
            iiStart=nXplore(kk);
            condLoop=true;
            pLoop=pLoop+1;
            while condLoop
                ii=nXplore(kk);
                nXplore(kk)=[];
                
                iip1=mod(ii,nX0)+1;
                if mod(ceil(iout(iip1)),nP1)~=mod(ceil(iout(ii)),nP1)
                    actTInd=ceil(mod(iout(ii),nP1));
                    actTInd(actTInd==0)=iout(ii);
                    [branchFollowProf,on]=inpolygon(profileCoord(actTInd,1),profileCoord(actTInd,2),...
                        targCoord(:,1),targCoord(:,2));
                else
                    actTInd=ceil(mod(jout(ii),nP2));
                    actTInd(actTInd==0)=jout(ii);
                    [branchFollowProf,on]=inpolygon(targCoord(actTInd,1),targCoord(actTInd,2),...
                        profileCoord(:,1),profileCoord(:,2));
                    branchFollowProf=~branchFollowProf;
                end
                
                if branchFollowProf
                    [ind1,iip1,branchCoord]=FollowBranchLoop(ii,iout,profileCoord,interCoord,nP1);
                else 
                    [ind1,iip1,branchCoord]=FollowBranchLoop(ii,jout,targCoord,interCoord,nP2);
                end
                newloop(pLoop).coord=[newloop(pLoop).coord;
                    branchCoord];
                kk=find(nXplore==iip1);
                condLoop=iip1~=iiStart && ~isempty(kk);
            end
           
        end
        newloop=newloop(1:pLoop);
    else
        inTarg=inpolygon(profileCoord(1,1),profileCoord(1,2),...
            targCoord(:,1),targCoord(:,2));
        
        inProf=inpolygon(targCoord(1,1),targCoord(1,2),...
            profileCoord(:,1),profileCoord(:,2));
        
        if ~inProf && ~inTarg
            newloop=[];
        elseif inProf && ~inTarg
            newloop.coord=targCoord;
        elseif ~inProf && inTarg
            newloop.coord=profCoord;
            
        else
            error('Intersections not detected but both profiles in the other')
        end
        
    end
    
end

function [ind1,iip1,branchCoord]=FollowBranchLoop(ii,iout,coord,interCoord,nP1)
    % Identifies the next intersection and returns the indices and
    % coordinates that correspond to the branch.
       
    
    [~,iip1]=min(mod(iout([1:ii-1,ii+1:end])-iout(ii),max(iout)));
    iip1=iip1+(iip1>=ii);
    
    
    if floor(iout(ii))~=floor(iout(iip1))
        if ceil(iout(ii))>floor(iout(iip1))
            ind1=[ceil(iout(ii)):nP1,1:floor(iout(iip1))];
        else
            ind1=ceil(iout(ii)):floor(iout(iip1));
        end
    else
        ind1=[];
    end
    ind1=mod(ind1-1,nP1)+1;
    branchCoord=[interCoord(ii,:);coord(ind1,:);interCoord(iip1,:)];
    
end

function [newloop]=XORProfileOP(profileCoord,targCoord)
    % compares two loops and returns the area not filled by both
    if ~CCWLoop(profileCoord)
        profileCoord=flip(profileCoord);
    end
    if ~CCWLoop(targCoord)
        targCoord=flip(targCoord);
    end
    
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
    
    nX0=numel(x0);
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    
    if nX0>0
        newloop=repmat(struct('coord',zeros(0,2)),[1,nX0]);
        nXplore=1:nX0;
        pLoop=0;
        while ~isempty(nXplore)
            % Need to make sure both profiles go in the same direction
            
            % THis code is insufficient it does not check for that no other
            % intersection lies on ind2 (cannot lie on ind1 as it is sorted)
            kk=1;
            iiStart=nXplore(kk);
            condLoop=true;
             pLoop=pLoop+1;
            while condLoop
                ii=nXplore(kk);
                nXplore(kk)=[];
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
                iij=ii;
                if floor(jout(ii))~=floor(jout(iip1))
                    if ceil(jout(ii))>floor(jout(iip1))
                        
                        isInTheWay1=(0<jout & jout<jout(iip1));
                        isInTheWay2=(jout>jout(ii) & jout<nP2);
                        if any(isInTheWay1)
                            subIsInTheWay=find(isInTheWay1);
                            [~,posInWay]=max(jout(subIsInTheWay));
                            iij=subIsInTheWay(posInWay);
                            ind2=flip(ceil(jout(iij)):floor(jout(iip1)));
                        elseif any(isInTheWay2)
                            subIsInTheWay=find(isInTheWay2);
                            [~,posInWay]=max(jout(subIsInTheWay));
                            iij=subIsInTheWay(posInWay);
                            ind2=flip([ceil(jout(iij)):nP2,1:floor(jout(iip1))]);
                        else
                            ind2=flip([ceil(jout(iij)):nP2,1:floor(jout(iip1))]); % iip1 does not change
                        end
                    else
                        
                        isInTheWay=jout>jout(ii) & jout<jout(iip1);
                        if any(isInTheWay)
                            subIsInTheWay=find(isInTheWay);
                            [~,posInWay]=max(jout(subIsInTheWay));
                            iij=subIsInTheWay(posInWay);
                        end
                        ind2=flip(ceil(jout(iij)):floor(jout(iip1)));
                        
                    end
                else
                    ind2=[];
                end
                
                
                
                newloop(pLoop).coord=[newloop(pLoop).coord;
                    [x0(ii),y0(ii)];
                    profileCoord(ind1,:);
                    [x0(iip1),y0(iip1)];
                    targCoord(ind2,:)];
                kk=find(nXplore==iij);
                condLoop=iij~=iiStart && ~isempty(kk);
            end
           
            
        end
        newloop=newloop(1:pLoop);
    else
        
        error('No Intersection, should not happen')
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end