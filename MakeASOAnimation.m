function [loopDatas]=MakeASOAnimation(surface,dataLog,initialLoops,targetLoops)
    
    
    plotPoints= @(points,c) plot(points([1:end,1],1),points([1:end,1],2),c);
    h=figure;
    ax=axes;
    
    loopDatas.static=struct('inital',initialLoops,'target',targetLoops);
    loopDatas.dynamic=struct('stage',[],'major',[],'surf',[],'ctrl',[]);
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
                    ll=ll+1;
                end
            end
        end
    end
    
    
    
end

function [cellTec]=TecplotOutputASOAnimation(loopDatas)
    cellTec=cell(0);
    
    fields=fieldnames(loopDatas.static);
    for ii=1:numel(fields)
        cellTec=[cellTec,CellEdgeMesh(loopDatas.static.(fields{ii}),[])];
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
    
    [loopsCtrl]=ControlPointLoops(ctrlPts,surface.Pn{nLevel});
    
    
    
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

function [loops]=ControlPointLoops(pts,Pn)
    [faces]=ConnectCtrlPts(Pn);
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