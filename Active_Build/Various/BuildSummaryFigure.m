



function []=BuildSummaryFigure(vecFig,cellNames)
    
    hSum=figure;
    axSum=axes;
    
    
    for ii=1:length(vecFig)
        ExtractFullLines(vecFig(ii),cellNames{ii},axSum)
    end
    UpdateLineColors(axSum)
end

function []=ExtractFullLines(figH,name,axSum)
    
    axh=findobj(figH,'type','axes');
    
    for ii=1:length(axh)
        if sum(axh(ii).View==[0 90])==2
            linesh=findobj(axh(ii),'type','line');
            kk=1;
            for jj=1:length(linesh)
                if strcmp(linesh(jj).LineStyle,'-')
                    newL=copyobj(linesh(jj),axSum);
                    newL.DisplayName=[name,' - ',newL.DisplayName,' - ', int2str(kk)];
                    kk=kk+1;
                end
            end
        end
    end
    
end

function []=UpdateLineColors(axSum)
    
    
    linesh=findobj(axSum,'type','line');
    cOrd=get(gca,'ColorOrder');
    for ii=1:length(linesh)
        linesh(ii).Color=cOrd(mod(ii-1,7)+1,:);
    end
    legend(linesh);
    
end

