function [h]=PlotConvMLSO(stageData)
    isLoadColor = false;
    h = figure;
    hold on
    [resstruct]=StageData2ResStruct(stageData);
    for ii=1:numel(resstruct.refine)
                    
        figobj.l(ii,1)=plot(resstruct.refine(ii).iterstart:resstruct.refine(ii).iterend,...
            resstruct.refine(ii).obj);
        
        if isLoadColor
            figobj.l(ii,1).Color=figobj.color;
        end
        figobj.l(ii,1).DisplayName=['Refine Stage ',int2str(resstruct.refine(ii).stage)];
        
        figobj.l(ii,2)=plot(resstruct.refine(ii).iterend,...
            resstruct.refine(ii).obj(end),'*','Color',figobj.l(ii,1).Color);
        figobj.l(ii,3)=plot([1,resstruct.refine(end).iterend+5],...
            resstruct.refine(ii).obj(end)*[1 1],'--','Color',figobj.l(ii,1).Color);
        
        figobj.t(ii,1)=text(resstruct.refine(end).iterend+8,...
            resstruct.refine(ii).obj(end),int2str(resstruct.refine(ii).ndesvar));
        figobj.t(ii,1).VerticalAlignment='middle';
        %figobj.t(ii,1).BackgroundColor=[1 1 1];
        figobj.t(ii,1).Color=figobj.l(ii,1).Color;
        
    end
end

function [resstruct]=StageData2ResStruct(stageData)
    kk = 1;
    totIter = sum([stageData.nMajIter]);
    nRef = numel(stageData);
    resstruct.refine = repmat(...
        struct('stage',[], 'iterstart', [], 'iterend', [], 'ndesvar',[],...
            'obj', [])...
        ,[1,nRef]);

    for ii = 1:nRef
        resstruct.refine(ii).stage =ii;
        resstruct.refine(ii).iterstart =kk;
        resstruct.refine(ii).iterend =kk+stageData(ii).nMajIter-1;
        kk = kk+stageData(ii).nMajIter;
        resstruct.refine(ii).ndesvar = numel(stageData(ii).majorData(1).x);
        resstruct.refine(ii).obj =[stageData(ii).majorData(:).objective];
        if ii>1
            resstruct.refine(ii).iterstart =resstruct.refine(ii).iterstart-1;
            resstruct.refine(ii).obj = [resstruct.refine(ii-1).obj(end),resstruct.refine(ii).obj];
        end
    end
end