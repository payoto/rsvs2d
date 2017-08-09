function []=ExploreFunction(func,numarg,cellarg,plotFuncs,nameCol,nameCell)
    
    if nargin<6
        nameCell=cellstr(num2str([1:size(cellarg,1)]'))';
        if nargin<5
            nameCol=0;
        end
    end
    lStyle={'-','--','-.'};
    if nameCol~=0
        numberCol=cellstr(num2str(numarg(:,nameCol),'%.1f'))';
    else
        numberCol=cellstr(num2str([1:size(numarg,1)]'))';
    end
    
    figure, 
    figC=ceil(sqrt(numel(plotFuncs)));
    figR=ceil(numel(plotFuncs)/figC);
    for ii=1:numel(plotFuncs)
        subplot(figR,figC,ii)
        hold on
    end
    
    numargcell=num2cell(numarg);
    numarg(:,max(nameCol,1))=numarg(:,max(nameCol,1))+1e-5;
    numargcellDel=num2cell(numarg);
    
    for ii=1:size(numarg,1)
        for jj=1:size(cellarg,1)
            pts=func(numargcell{ii,:},cellarg{jj,:})';
            ptsDel=func(numargcellDel{ii,:},cellarg{jj,:})';
            [~,ptsNorm]=NormalDistance(pts,ptsDel);
            nameLine=[nameCell{jj},':',numberCol{ii}];
            for kk=1:numel(plotFuncs)
                subplot(figR,figC,kk)
                
                l{kk}(ii,jj)=plotFuncs{kk}(pts,ptsDel,ptsNorm);
                
                l{kk}(ii,jj).DisplayName=nameLine;
                l{kk}(ii,jj).Color=l{kk}(ii,1).Color;
                l{kk}(ii,jj).LineStyle=lStyle{mod(jj-1,numel(lStyle))+1};
            end
        end
    end
    
    for ii=1:numel(plotFuncs)
        legend(l{ii}(:))
    end
    
    
end