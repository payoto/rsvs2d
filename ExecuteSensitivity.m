

%% Sensitivity method development

function [snaxmode]=ExecuteSensitivity(snaxel,snakposition,sensSnax,volumefraction,isPlot)
    
    [snaxmode]=ExtractSensitivity(snaxel,snakposition,sensSnax,volumefraction,isPlot);
    
    rootloop=OrderSurfaceSnaxel(snaxel);
    for jj=1:numel(rootloop)
        A(jj)=abs(CalculatePolyArea(rootloop(jj).snaxel.coord));
    end
    startVol=sum(A);
    
    for ii=1:numel(snaxmode)
        for jj=1:numel(snaxmode(ii).loopsnaxel)
            snaxmode(ii).A(jj)=abs(CalculatePolyArea(snaxmode(ii).loopsnaxel(jj).snaxel.coord));
        end
        [snaxmode(ii).vecmode]=BuildVectors(rootloop,snaxmode(ii).loopsnaxel);
        snaxmode(ii).deltaFrac=(sum(snaxmode(ii).A)-startVol)/snaxmode(ii).cellVol;
    end
end

function [vectormode]=BuildVectors(looproot,loopsnaxel)
    
    for ii=1:numel(loopsnaxel)
        actRoot=FindObjNum([],loopsnaxel(ii).snaxel.index,looproot(ii).snaxel.index);
        subAct=SubDivision(loopsnaxel(ii).snaxel.coord,2,'chaikin',[0 0],'none');
        subRoot{ii}=SubDivision(looproot(ii).snaxel.coord(actRoot,:),2,'chaikin',[0 0],'none');
        vecMode{ii}=subAct-subRoot{ii};
    end
    vectormode.points=vertcat(subRoot{:});
    vectormode.vector=vertcat(vecMode{:});
end

function [snaxmode]=ExtractSensitivity(snaxel,snakposition,sensSnax,volumefraction,isPlot)
    
    if nargin<5
       isPlot=false; 
    end
    global unstructglobal
    t1=now;
    cellVols=[volumefraction(find(sum(abs(sensSnax))~=0)).totalvolume];
    cellInd=[volumefraction(find(sum(abs(sensSnax))~=0)).oldCellInd];
    sensSnax(:,find(sum(abs(sensSnax))==0))=[];
    dCurr=[snaxel(:).d];
    
    [snaxOrd]=SplitSnaxLoops(snaxel); % Isolates individual loops
    maxDistRatio=1/10;
    [dChange]=FindModalDistanceChange(sensSnax,maxDistRatio);
    
    
    coord1=vertcat(snakposition(:).coord);
    [dir]=sum((vertcat(snakposition(:).vector)~=0).*[ones([length(snaxel), 1]),...
        ones([length(snaxel), 1])*2],2);
    
    for ii=1:length(sensSnax(1,:)),
        snaxCopy(ii,:)=snaxel;
        
        dAct=dCurr'+dChange{ii};
        for jj=1:length(snaxel)
            snaxCopy(ii,jj).d=dAct(jj);
        end
        [snakposition2(ii,:)]=PositionSnakes(snaxCopy(ii,:),unstructglobal);
    end
    
    coordBase{numel(snaxOrd)}=[];
    refPos{numel(snaxOrd)}=[];
    coordNew{numel(snaxOrd),numel(dChange)}=[];
    newPos{numel(snaxOrd),numel(dChange)}=[];
    newOrd{numel(snaxOrd),numel(dChange)}=[];
    for ii=1:numel(snaxOrd)
        coordBase{ii}=vertcat(snakposition(snaxOrd{ii}).coord);
        [refPos{ii}]=CreateComparisonMatrix(coordBase{ii});
        for jj=1:numel(dChange)
            coordNew{ii,jj}=vertcat(snakposition2(jj,snaxOrd{ii}).coord);
            [newPos{ii,jj}]=CreateComparisonMatrix(coordNew{ii,jj});
            [newOrd{ii,jj}]=(CompareTestpos(refPos{ii},newPos{ii,jj},1:numel(snaxOrd{ii}),dir(snaxOrd{ii})));
        end
    end
    
    for jj=1:numel(dChange)
        snaxmode(jj).cellindex=cellInd(jj);
        snaxmode(jj).snaxel=snaxCopy(jj,:);
        snaxOrdTemp=newOrd(:,jj);
        for ii=1:numel(snaxOrd)
            snaxOrdTemp{ii}=snaxOrd{ii}(snaxOrdTemp{ii});
        end
        [snaxmode(jj).snaxel]=MatchOrders(snaxmode(jj).snaxel,snaxOrdTemp);
        snaxmode(jj).loopsnaxel=OrderSurfaceSnaxel(snaxmode(jj).snaxel);
        snaxmode(jj).cellVol=cellVols(jj);
    end
    
    
    t2=now;
    
    
    if isPlot
        for jj=1:numel(dChange)
            coordFull{jj}=vertcat(snakposition2(jj,:).coord);
        end
        
        for jj=1:numel(dChange)
            figure
            hold on
            for ii=1:numel(snaxOrd)
                
                plot(coordBase{ii}(:,1),coordBase{ii}(:,2),'+-')
                plot(coordNew{ii,jj}(newOrd{ii,jj},1),coordNew{ii,jj}(newOrd{ii,jj},2),'o-')
                
            end
            for ii=1:size(coord1,1)
                plot([coord1(ii,1),coordFull{jj}(ii,1)],[coord1(ii,2),coordFull{jj}(ii,2)],'k--')
            end
            
        end
    end
    t3=now;
    %         coord2=vertcat(snakposition2(:).coord);
    %         [testPos2]=CreateComparisonMatrix(coord2);
    %         [newOrd]=CompareTestpos(testPos1,testPos2,snaxOrd,dir);
    %
    
    %         figure
    %         plot(coord1(snaxOrd,1),coord1(snaxOrd,2),'o-',coord2(newOrd,1),coord2(newOrd,2),'o-')
    %         hold on
    %         for jj=1:length(newOrd)
    %             plot([coord1(newOrd(jj),1),coord2(newOrd(jj),1)],[coord1(newOrd(jj),2),coord2(newOrd(jj),2)],'k--')
    %         end
    %         title(['mode ',int2str(ii)])
    %         axis equal
    
    
    datestr(t3-t2,'SS:FFF')
    datestr(t2-t1,'SS:FFF')
    
end

function [snaxel]=MatchOrders(snaxel,snaxOrd)
    
    for ii=1:numel(snaxOrd)
        
        nOrd=numel(snaxOrd{ii});
        for jj=1:nOrd
            ind=mod(jj-1,nOrd)+1;
            indp1=mod(jj,nOrd)+1;
            indm1=mod(jj-2,nOrd)+1;
            
            snaxel(snaxOrd{ii}(ind)).snaxprec=snaxel(snaxOrd{ii}(indm1)).index;
            snaxel(snaxOrd{ii}(ind)).snaxnext=snaxel(snaxOrd{ii}(indp1)).index;
            snaxel(snaxOrd{ii}(ind)).connectivity=...
                [snaxel(snaxOrd{ii}(indm1)).index,snaxel(snaxOrd{ii}(indp1)).index];
        end
        
    end
    
    snaxel=snaxel([snaxOrd{:}]);
    
end

function [snakposition]=PositionSnakes(snaxel,unstructured)
    % Returns an array with Snaxel coordinates preceded by snaxel indices
    vertIndex=unstructured.vertex.index;
    vertCoord=unstructured.vertex.coord;
    fromVertex=[snaxel(:).fromvertex];
    toVertex=[snaxel(:).tovertex];
    
    nSnaxel=length(snaxel);
    
    for ii=nSnaxel:-1:1
        iToVert=vertCoord(find(vertIndex==toVertex(ii)),:); %#ok<FNDSB> % extract vertex coordinates
        iFromVert=vertCoord(find(vertIndex==fromVertex(ii)),:); %#ok<FNDSB>
        
        snakposition(ii).index=snaxel(ii).index;
        snakposition(ii).coord=iFromVert+(iToVert-iFromVert)*snaxel(ii).d;
        snakposition(ii).vectornotnorm=(iToVert-iFromVert);
        snakposition(ii).vertInit=iFromVert;
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

function [snaxOrd]=SplitSnaxLoops(snaxel)
    % Splits a snake into its component loops
    
    kk=1;
    jj=1;
    snaxInd=[snaxel(:).index];
    %snaxOrd=zeros(size(snaxel));
    ordList=1:numel(snaxel);
    ll=1;
    for ii=1:length(snaxel)
        kk=FindObjNum([],[snaxel(kk).snaxnext],snaxInd);
        if ordList(kk)==0
            jj=jj+1;
            kk=min(ordList(ordList~=0));
            ll=1;
        end
        snaxOrd{jj}(ll)=kk;
        ordList(kk)=0;
        ll=ll+1;
    end
    
end

function [dChange]=FindModalDistanceChange(sensSnax,maxDistRatio)
    
    nMode=size(sensSnax,2);
    dChange{nMode}=[];
    for ii=1:nMode
        dChange{ii}=sensSnax(:,ii)/max(abs(sensSnax(:,ii)))*maxDistRatio;
    end
end

function []=testSensitivity(snaxel,snakposition,sensSnax)
    
    global unstructglobal
    sensSnax(:,find(sum(abs(sensSnax))==0))=[];
    dCurr=[snaxel(:).d];
    kk=1;
    snaxInd=[snaxel(:).index];
    snaxOrd=zeros(size(snaxel));
    ordList=1:numel(snaxOrd);
    for ii=1:length(snaxel)
        kk=FindObjNum([],[snaxel(kk).snaxnext],snaxInd);
        if ordList(kk)==0
            kk=min(ordList(ordList~=0));
        end
        snaxOrd(ii)=kk;
        ordList(kk)=0;
    end
    %snaxOrd(end+1)=snaxOrd(1);
    coord1=vertcat(snakposition(:).coord);
    [dir]=sum((vertcat(snakposition(:).vector)~=0).*[ones([length(snaxel), 1]),...
        ones([length(snaxel), 1])*2],2);
    [testPos1]=CreateComparisonMatrix(coord1);
    l=max(sum(vertcat(snakposition(:).vectornotnorm).^2,2));
    for ii=1:length(sensSnax(1,:)),
        snaxCopy=snaxel;
        e1=(1)./sensSnax(:,ii);
        e2=(-1)./sensSnax(:,ii);
        e1_sel=min(e1(e1>0));
        e2_sel=min(e2(e2>0));
        e_sel(ii)=min([e1_sel,e2_sel]);
        dChange{ii}=sensSnax(:,ii)/max(abs(sensSnax(:,ii)))*1;
        %dChange{ii}=-sensSnax(:,ii)/100;
        dAct=dCurr'+dChange{ii};
        
        for jj=1:length(snaxel)
            snaxCopy(jj).d=dAct(jj);
        end
        [snakposition2]=PositionSnakes(snaxCopy,unstructglobal);
        
        
        coord2=vertcat(snakposition2(:).coord);
        Delta{ii}=coord2-coord1;
        [testPos2]=CreateComparisonMatrix(coord2);
        [newOrd]=CompareTestpos(testPos1,testPos2,snaxOrd,dir);
        %newOrd=snaxOrd;
        figure
        plot(coord1(snaxOrd,1),coord1(snaxOrd,2),'+-',coord2(newOrd,1),coord2(newOrd,2),'o-')
        hold on
        for jj=1:length(newOrd)
            plot([coord1(newOrd(jj),1),coord2(newOrd(jj),1)],[coord1(newOrd(jj),2),coord2(newOrd(jj),2)],'k--')
        end
        title(['mode ',int2str(ii)])
        
    end
    
    
end

function [testPos]=CreateComparisonMatrix(coord)
    
    [m,n]=size(coord);
    testPos{n}=[];
    for jj=1:n
        testPos{jj}=zeros(m);
        for ii=1:m
            testPos{jj}(:,ii)=(coord(ii,jj)>coord(:,jj))+(coord(ii,jj)<coord(:,jj))*-1;
        end
    end
    
    
end

function [newOrd]=CompareTestpos(t1,t2,ord,dir)
    
    oldOrd=ord;
    ord=FindObjNum([],ord,ord)';
    
    for ii=1:length(t1)
        tD{ii}=t1{ii}~=t2{ii};
    end
    tX=true(size(tD{1}));
    tDel=zeros(size(tD{1}));
    for ii=1:length(tD)
        tX=tX & tD{ii};
        
    end
    for ii=1:length(tD)
        tDel=tDel+ii*(tD{ii} & ~tX);
    end
    
    [newOrd,tX]=ReorderList(ord,tX);
    newOrd=newOrd(end:-1:1);
    [newOrd,tX]=ReorderList(newOrd,tX);
    newOrd=newOrd(end:-1:1);
    
    %newOrd=[newOrd(end-1),newOrd];
    rmSnak=[];
    kk=1;
    for ii=1:length(newOrd)
        
        ind=mod(ii-1,length(newOrd))+1;
        indp1=mod(ii,length(newOrd))+1;
        indm1=mod(ii-2,length(newOrd))+1;
        
        
        tTest=tDel([newOrd(indm1),newOrd(ind),newOrd(indp1)],[newOrd(ind)]);
        tTest2=t2{abs(dir(newOrd(ind))-3)}([newOrd(indm1),newOrd(ind),newOrd(indp1)]...
            ,[newOrd(ind)]);
        
        if (sum(tTest~=dir(newOrd(ind)) & tTest~=0)) ... % overtaken by neighbour
                || (sum(tTest2==0 & (tTest==dir(newOrd(ind))))) % neighbours on same line have crossed
            rmSnak(kk)=ind;
            kk=kk+1;
            
            
        end
    end
    newOrd(rmSnak)=[];
    newOrd=oldOrd(newOrd);
end

function [newOrd,tX]=ReorderList(ord,tX)
    
    newOrd=ord;
    for ii=0:2*length(ord)-1
        ind1=mod(ii,length(ord))+1;
        ind2=mod(ii+1,length(ord))+1;
        if tX(newOrd(ind2),newOrd(ind1)) % if crossing
            
            tX(newOrd(ind2),newOrd(ind1))=false; % delete the crossing flag
            tX(newOrd(ind1),newOrd(ind2))=false;
            
            interim=newOrd(ind1); % invert connection
            newOrd(ind1)=newOrd(ind2);
            newOrd(ind2)=interim;
        end
    end
    
    
end

function [loopsnaxel]=OrderSurfaceSnaxel(snaxel)
    % function extracting the snaxels into their separate loops
    global unstructglobal
    
    snaxPositions=PositionSnakes(snaxel,unstructglobal);
    
    nSnax=length(snaxel);
    blockSegments=zeros(2*nSnax,2);
    for ii=1:nSnax
        for jj=0:1
            blockSegments(2*ii-jj,:)=[snaxel(ii).index,snaxel(ii).connectivity(jj+1)];
        end
    end
    
    cellSimilar=FindIdenticalVector(blockSegments);
    for ii=1:length(cellSimilar)
        blockEdgeIndex(ii)=cellSimilar{ii}(1);
    end
    blockEdges=blockSegments(blockEdgeIndex,:);
    % Order edges into closed loops
    [cellOrderedVertex]=OrderBlockEdges(blockEdges);
    snaxIndex=[snaxel(:).index];
    for ii=1:length(cellOrderedVertex)
        loopsnaxel(ii).snaxel.index=[cellOrderedVertex{ii}(:,1)];
        loopIndices=FindObjNum(snaxel,loopsnaxel(ii).snaxel.index,snaxIndex);
        loopsnaxel(ii).snaxel.coord=vertcat(snaxPositions(loopIndices).coord);
        %loopsnaxel(ii).edge.index=isEdgeIndex(cellOrderedEdges{ii});
    end
    
end

function [A]=CalculatePolyArea(points)
    
    pointsVec=points';
    pointsVec=pointsVec(:);
    %plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;
    
    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
    
end
