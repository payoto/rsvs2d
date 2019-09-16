
function [A, newloop]=IntersectingArea(profileCoord,targCoord)
    % compares two loops and returns the area filled by both
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
        keepLoop = true(size(newloop));
        
        for ii=1:pLoop
            keepLoop = size(newloop(ii).coord,1)>2;
        end
        newloop=newloop(keepLoop);
        
        
        
    else
        % Is one in the other?
        newloop(1).coord = targCoord;
        newloop(2).coord = profileCoord;
        
        [in1,on1] = inpolygon(profileCoord(:,1),profileCoord(:,2),...
            targCoord(:,1),targCoord(:,2));
        if any(in1) && ~all(in1)
            error('Contradiction')
        end
        if any(in1)
            newloop(1).coord = targCoord;
            newloop = newloop(1);
        end
        
        if ~any(in1)
            [in1,on1] = inpolygon(targCoord(:,1),targCoord(:,2),...
                profileCoord(:,1),profileCoord(:,2));
            if any(in1) && ~all(in1)
                error('Contradiction')
            end
            if any(in1)
                newloop(1).coord = profileCoord;
                newloop = newloop(1);
            end
        end
        
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end