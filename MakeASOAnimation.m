function [loopDatas,cellTec]=MakeASOAnimation(surface,dataLog,initialLoops,targetLoops)
    
    
    plotPoints= @(points,c) plot(points([1:end,1],1),points([1:end,1],2),c);
    h=figure;
    ax=axes;
    cellTec=cell(0);
    
    if ~isstruct(surface)

        loopDatas.static=struct('initial',[],'target',[],'objective',[]);
        loopDatas.static.initial=initialLoops;
        loopDatas.static.target=targetLoops;
        loopDatas.dynamic=struct('stage',[],'major',[],'surf',[],'ctrl',[],'obj',[]);
        ll=1;
        [surfLoop]=OrderBlockEdges(surface.geom.faces);
        for ii=1:numel(dataLog.stageData)
            if ~isempty(dataLog.stageData(ii).nMajIter)
                for jj=1:numel(dataLog.stageData(ii).majorData)
                    if ~isempty(dataLog.stageData(ii).majorData(jj).xStep)
                        cla(ax,'reset')
                        hold(ax,'on')
                        title(['Level ',int2str(ii),', Major ',int2str(jj)])
                        for kk=1:numel(initialLoops)
                            plotPoints(initialLoops{kk},'-');
                            ax.ColorOrderIndex=ax.ColorOrderIndex-1;
                        end
                        ax.ColorOrderIndex=ax.ColorOrderIndex+1;
                        for kk=1:numel(targetLoops)
                            plotPoints(targetLoops{kk},'-');
                            ax.ColorOrderIndex=ax.ColorOrderIndex-1;
                        end
                        ax.ColorOrderIndex=ax.ColorOrderIndex+1;
                        [loopDatas.dynamic(ll).stage,loopDatas.dynamic(ll).surf,...
                            loopDatas.dynamic(ll).ctrl]=PlotSingleStep(ax,...
                            surface,dataLog.stageData(ii).majorData(jj),surfLoop);
                        loopDatas.dynamic(ll).major=jj;
                        loopDatas.static.objective(ll)=dataLog.stageData(ii).majorData(jj).objective;
                        loopDatas.dynamic(ll).obj=dataLog.stageData(ii).majorData(jj).objective;
                        ll=ll+1;
                    end
                end
            end
        end
        
    else
        loopDatas=surface;
    end
    assignin('base','loopDatas',loopDatas)
    [cellTec]=TecplotOutputASOAnimation(loopDatas);
end

function [cellTec]=TecplotOutputASOAnimation(loopDatas)
    cellTec=cell(0);
    n=ceil(log10(numel(loopDatas.dynamic)));
    for jj=1:numel(loopDatas.dynamic)
        t(jj)=loopDatas.dynamic(jj).stage+loopDatas.dynamic(jj).major/10^n;
        n2=sum([loopDatas.dynamic.stage]==loopDatas.dynamic(jj).stage)+1;
        t2(jj)=loopDatas.dynamic(jj).stage+loopDatas.dynamic(jj).major/n2;
    end
   
    cellMesh=OutputLineAsFelineSeg(loopDatas.static.objective,t,t2,loopDatas.dynamic);
    
    fields={'initial','target'};
    for ii=1:numel(fields)
        cellTec=[cellTec,Loop2CellEdgeMesh(loopDatas.static.(fields{ii}),[])];
    end
    
    fields={'surf','ctrl'};
    n=ceil(log10(numel(loopDatas.dynamic)));
    for jj=1:numel(loopDatas.dynamic)
        t=loopDatas.dynamic(jj).stage+loopDatas.dynamic(jj).major/10^n;
        for ii=1:numel(fields)
            cellTec=[cellTec,Loop2CellEdgeMesh(loopDatas.dynamic(jj).(fields{ii}),[],ii+5,t)];
        end
    end
    
    
    cellTec=[cellTec,cellMesh];
    
end

function cellMesh=OutputLineAsFelineSeg(objective,t,t2,dataDynamic)
    
    cellMesh=cell(0);
    vertDat=[t;objective]';
    velDat=[log10(objective);1:numel(t);t2]';
    connDat=[1:numel(objective)-1;2:numel(objective)]';
    [cellMesh]=CellEdgeMesh(vertDat,[1:numel(objective)],connDat,velDat);
    
    for ii=1:numel(dataDynamic)
        [cellMesh]=[cellMesh,CellEdgeMesh([t(ii),dataDynamic(ii).obj],...
            [1],[1 1],[log10(dataDynamic(ii).obj) ii t2(ii)],8,t(ii))];
    end
    
end

function [nLevel,loopSurf,loopsCtrl]=PlotSingleStep(ax,surface,dataStep,surfLoop)
    
    hold(ax,'on')
    ax.XLim=[-0.8 1.2];
    ax.YLim=[-0.13 0.13];
    plotPoints= @(points,c) plot(points([1:end,1],1),points([1:end,1],2),c);
    for ii=1:numel(surfLoop)
        loopSurf{ii}=dataStep.surface(surfLoop{ii}(:,1),1:2);
    end
    
    ctrlPts=[dataStep.x(1:(numel(dataStep.x)*0.5-2)),...
        dataStep.x(round(numel(dataStep.x)*0.5)+1:end-2)];
    nCtrl=size(ctrlPts,1);
    nLevel=find(cellfun(@(x)size(x,2),surface.phi_n)==nCtrl);
    
    [loopsCtrl]=ControlPointLoops(ctrlPts,surface.Pn{nLevel},surface.Pn,nLevel);
    
    
    
    for ii=1:numel(loopSurf)
        plotPoints(loopSurf{ii},'+-');
        ax.ColorOrderIndex=ax.ColorOrderIndex-1;
    end
    ax.ColorOrderIndex=ax.ColorOrderIndex+1;
    for ii=1:numel(loopsCtrl)
        plotPoints(loopsCtrl{ii},'s-');
        ax.ColorOrderIndex=ax.ColorOrderIndex-1;
    end
    ax.ColorOrderIndex=ax.ColorOrderIndex+1;
    axis equal
    pause(0.001)
end

function [loops]=ControlPointLoops(pts,Pn,PnArray,nLevel)
    [faces]=ConnectCtrlPts(Pn);
    if size(faces,1)~=size(Pn,2)
        [faces]=ConnectCtrlPtsFallBack(PnArray,nLevel);
    end
    [cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(faces);
    for ii=1:numel(cellOrderedVertex)
        loops{ii}=pts(cellOrderedVertex{ii}(:,1),1:2);
    end
end

function [faces]=ConnectCtrlPts(Pn)
    Pn=Pn>0;
    Pn=Pn(find(sum(Pn,2)==2),:);
    [i,j]=find(Pn>0);
    
    [~,iOrd]=sort(i);
    faces=j(iOrd);
    
    faces=reshape(faces,[2,numel(faces)/2])';
    
    
    
end


function [faces]=ConnectCtrlPtsFallBack(PnArray,nLevel)
  
   [loopInd]=ConnectCtrlPtsFallBackLink(PnArray{1});
   
   for ii=1:nLevel-1
       loopTransfer=zeros(size(PnArray{ii},1),1);
        for jj=1:size(PnArray{ii},1)
            loopTransfer(jj)=find(PnArray{ii}(jj,:),1,'first');
        end
        loopInd=loopInd(loopTransfer);
   end
            
    faces=zeros(0,2);
    
    for ii=1:max(loopInd);
        currInd=sort(find(loopInd==ii));
        faces=[faces;[currInd,currInd([2:end,1])]];
    end
    
end


function [loopInd]=ConnectCtrlPtsFallBackLink(Pn)
    Pn=Pn>0;
    
    loopInd=zeros(size(Pn,2),1);
    currLoop=1;    
    currInd=1;
    while any(loopInd==0)
        loopInd(currInd)=currLoop;
        [i,~]=find(Pn(:,currInd));
        Pn(:,currInd)=0;
        [~,currInd]=find(Pn(i,:));
        
        if isempty(currInd)
            currInd=find(loopInd==0,1,'first');
            currLoop=currLoop+1;
        end
        
    end
%     faces=zeros(0,2);
%     
%     for ii=1:currLoop-1;
%         currInd=sort(find(loopInd==ii));
%         faces=[faces;[currInd,currInd([2:end,1])]];
%     end
    
    
    
end