function [snaxmove]=BuildMovementStructures(snaxel,snakposition,snaxmode,sensSnax,volumefraction)
    
    if numel(snaxel)~=numel(snakposition)
        error('Dimensiom mismatch between snakposition and snaxel while building movement structure')
    elseif any([snaxel(:).index]~=[snakposition(:).index])
        error('Indices mismatch between snakposition and snaxel while building movement structure')
    end
    
    [sensSnax]=ScaleTrimSensSnax(sensSnax,volumefraction,snaxmode);
    [loopsnaxel,nSnax,nLoop]=CCWLoopLength(snaxel,snakposition);
    
    [snaxmove]=BuildSnaxMoveTemplate(nSnax,nSnax,nSnax);
    for ii=1:nLoop
        [snaxmove(ii)]=CalculateMoveData(snaxmove(ii),snakposition,sensSnax,...
            loopsnaxel(ii).snaxel);
    end
    deltaFill=zeros([3,numel(snaxmode)]);
    deltaFill(1,3)=1;deltaFill(2,10)=1;deltaFill(3,25)=1;
    deltaFill=deltaFill*1e-2;
    MoveToFill(snaxmove,deltaFill)
end

function [loopsnaxel,nSnax,nLoop]=CCWLoopLength(snaxel,snakposition)
    
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel,snakposition);
    nLoop=numel(loopsnaxel);
    nSnax=zeros([1,nLoop]);
    for ii=1:nLoop
        if ~CCWLoop(loopsnaxel(ii).snaxel.coord);
            loopsnaxel(ii).snaxel.index=flip(loopsnaxel(ii).snaxel.index);
            loopsnaxel(ii).snaxel.coord=flip(loopsnaxel(ii).snaxel.coord);
        end
        % Add the length coordinates
        loopsnaxel(ii).snaxel.edgelength=sqrt(sum((loopsnaxel(ii).snaxel.coord([2:end,1],:)-...
            loopsnaxel(ii).snaxel.coord).^2,2));
        loopsnaxel(ii).snaxel.lpos=cumsum([0;loopsnaxel(ii).snaxel.edgelength]);
        nSnax(ii)=numel(loopsnaxel(ii).snaxel.index);
    end
end

function [snaxmove]=BuildSnaxMoveTemplate(nSnax,nEdge,nVertex)
    
    nLoop=numel(nSnax);
    
    snaxmove=struct('snax',struct([]),'vertex',struct([]),'edge',struct([]),...
        'support',struct('nSnax',[],'maxL',[]));
    snaxmove=repmat(snaxmove,[1,nLoop]);
    
    snaxTemp=struct('index',[],'rN',[],'rT',[],'d',[],'sens',[],'posL',[]);
    edgeTemp=struct('index',0,'posL',[],'coord',[],'normal',[]);
    vertTemp=edgeTemp;
    
    for ii=1:nLoop
        
        snaxmove(ii).snax=repmat(snaxTemp,[1,nSnax(ii)]);
        snaxmove(ii).vertex=repmat(edgeTemp,[1,nVertex(ii)]);
        snaxmove(ii).edge=repmat(vertTemp,[1,nEdge(ii)]);
        snaxmove(ii).support.nSnax=nSnax(ii);
    end
    
end

function [snaxmove]=CalculateMoveData(snaxmove,snakposition,sensSnax,loopsnaxel)
    
    fullSnax=[snakposition(:).index];
    
    snaxSub=FindObjNum([],loopsnaxel.index,fullSnax);
    snaxSub=reshape(snaxSub,[1,numel(snaxSub)]);
    
    snaxmove.support.maxL=loopsnaxel.lpos(end);
    
    
    kk=1;
    for ii=snaxSub
        kks=[kk:min([kk+1,snaxmove.support.nSnax]),...
            mod(kk,snaxmove.support.nSnax)+1:1];
        snaxmove.snax(kk).index=snakposition(ii).index;
        snaxmove.snax(kk).d=sqrt(sum(snakposition(ii).vectornotnorm.^2));
        snaxmove.snax(kk).sens=sensSnax(ii,:);
        snaxmove.snax(kk).posL=loopsnaxel.lpos(kk);
        
        snaxmove.vertex(kk).index=snakposition(ii).index;
        snaxmove.vertex(kk).normal=sum(snakposition(ii).normvector);
        snaxmove.vertex(kk).normal=snaxmove.vertex(kk).normal/...
            sqrt(sum(snaxmove.vertex(kk).normal.^2));
        snaxmove.vertex(kk).coord=snakposition(ii).coord;
        snaxmove.vertex(kk).posL=loopsnaxel.lpos(kk);
        
        snaxmove.edge(kk).index=loopsnaxel.index(kks);
        snaxmove.edge(kk).coord=mean(loopsnaxel.coord(kks,:));
        snaxmove.edge(kk).posL=mean(loopsnaxel.lpos(kk:kk+1));
        snaxmove.edge(kk).normal=snakposition(ii).vectornext;
        
        [vecAngles]=ExtractAnglepm180(snaxmove.vertex(kk).normal,snakposition(ii).vector);
        snaxmove.snax(kk).rT=sin(vecAngles);
        snaxmove.snax(kk).rN=cos(vecAngles);
        kk=kk+1;
    end
    
    
    
    
end

function [sensSnax]=ScaleTrimSensSnax(sensSnax,volumefraction,snaxmode)
    % Scale the sensitivity values in sensSnax to be in units of d per unit
    % of volume fraction
    
    
    cellInd=[volumefraction(:).oldCellInd];
    cellSens=[snaxmode(:).cellindex];
    
    sensRatio=[snaxmode(:).sensSnaxRatio]./[snaxmode(:).deltaFrac];
    colSub=FindObjNum([],cellSens,cellInd);
    
    sensSnax=sensSnax(:,colSub).*(ones([size(sensSnax,1),1])*sensRatio);
    
    
end

% Movement function

function [newloop]=MoveToFill(snaxmove,deltaFill)
    % deltaFill is a matrix where each row is a different fill
    
    
    % Find new Lpos
    
    nFill=size(deltaFill,1);
    newloop=repmat(struct('snaxel',struct('index',[],'coord',[])),[nFill,numel(snaxmove)]);
    
    % deltaFill must be trimmed to only have active snax boxes
    for ii=1:numel(snaxmove)
        % each column represent a different fill
        % each row represents a different snaxel
        
        % preparation
        nSnax=numel(snaxmove(ii).snax);
        nEdge=numel(snaxmove(ii).edge);
        nVert=numel(snaxmove(ii).vertex);
        snaxList=[snaxmove(ii).snax(:).index];
        vertPosL=repmat(reshape([snaxmove(ii).vertex(:).posL],...
            [1 1 nVert]),[nSnax,nFill,1]);
        edgePosL=repmat(reshape([snaxmove(ii).edge(:).posL],...
            [1 1 nEdge]),[nSnax,nFill,1]);
        
        % vertex recalc
        mvtCoeffFill=(vertcat(snaxmove(ii).snax(:).sens)*deltaFill')...
            .*(vertcat(snaxmove(ii).snax(:).d)*(ones([1,size(deltaFill,1)])));
        
        coeffNormpos=mvtCoeffFill.*(vertcat(snaxmove(ii).snax(:).rN)*(ones([1,size(deltaFill,1)])));
        
        deltaLpos=mvtCoeffFill.*(vertcat(snaxmove(ii).snax(:).rT)*(ones([1,size(deltaFill,1)])));
        newLpos=deltaLpos+(vertcat(snaxmove(ii).snax(:).posL)*(ones([1,size(deltaFill,1)])));
        newLpos(newLpos<0)=snaxmove(ii).support.maxL+newLpos(newLpos<0);
        newLpos(newLpos>snaxmove(ii).support.maxL)=...
            newLpos(newLpos<0)-snaxmove(ii).support.maxL;
        [~,newOrder]=sort(newLpos);
        
        [distVert,closeVert]=min(abs(vertPosL-repmat(newLpos,[1 1 nVert])),[],3);
        [distEdge,closeEdge]=min(abs(edgePosL-repmat(newLpos,[1 1 nEdge])),[],3);
        
        for jj=1:nFill
            
            vertCoord=reshape(vertcat(snaxmove(ii).vertex(closeVert(:,jj)).coord),[nVert,2,1]);
            edgeCoord=reshape(vertcat(snaxmove(ii).edge(closeEdge(:,jj)).coord),[nVert,2,1]);
            newCoordTemp=(vertCoord.*(distEdge(:,jj)*ones([1 2]))...
                +edgeCoord.*(distVert(:,jj)*ones([1 2])))...
                ./((distVert(:,jj)+distEdge(:,jj))*ones([1 2]));
            newLCoord=newCoordTemp;
            
            vertNormal=reshape(vertcat(snaxmove(ii).vertex(closeVert(:,jj)).normal),[nVert,2,1]);
            edgeNormal=reshape(vertcat(snaxmove(ii).edge(closeEdge(:,jj)).normal),[nVert,2,1]);
            newNormalTemp=(vertNormal.*(distEdge(:,jj)*ones([1 2]))...
                +edgeNormal.*(distVert(:,jj)*ones([1 2])))...
                ./((distVert(:,jj)+distEdge(:,jj))*ones([1 2]));
            newNormal=newNormalTemp./(sqrt(sum(newNormalTemp.^2,2))*ones([1 2]));
            
            newCoord=(coeffNormpos(:,jj)*ones([1 2])).*...
                newNormal+newLCoord;
            
            newloop(jj,ii).snaxel.index=snaxList(newOrder(:,jj));
            newloop(jj,ii).snaxel.coord=newCoord(newOrder(:,jj),:);
        end

    end
    
    
    if false
    
        plotPoints=@(points,str) plot(points([1:end,1],1),points([1:end,1],2),str);
        quiverPoints= @(points,vector,str) quiver(points([1:end,1],1),...
            points([1:end,1],2),vector([1:end,1],1),vector([1:end,1],2));
        figure,hold on
        axis equal


        colorStr='rky';

        for jj=1:nFill
            for ii=1:numel(snaxmove),
                plotPoints(newloop(jj,ii).snaxel.coord,[colorStr(jj),'-s'])
            end
        end
        for ii=1:numel(snaxmove),
            plotPoints(vertcat(snaxmove(ii).vertex(:).coord),'b-+')
            quiverPoints(vertcat(snaxmove(ii).vertex(:).coord),vertcat(snaxmove(ii).vertex(:).normal))
            plotPoints(vertcat(snaxmove(ii).edge(:).coord),'go')
            quiverPoints(vertcat(snaxmove(ii).edge(:).coord),vertcat(snaxmove(ii).edge(:).normal))
        end
    end
end