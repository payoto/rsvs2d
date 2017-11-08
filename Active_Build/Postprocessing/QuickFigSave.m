function []=QuickFigSave(figList)
    
    addpath ..\..\matlabpackages\export_fig
    for ii=figList;
        h=findobj(ii);
        h(1).Renderer='opengl';
        h(1).Color='none';
        axh=findobj(h(1),'type','axes');
        for jj=1:numel(axh)
            hCol{jj}=axh(jj).Color;
            axh(jj).Color='none';
        end
        figName=matlab.lang.makeValidName([h(1).Name,int2str(ii)]);
        %print(h(1),'-r300','-dpng',['.\fig\',figName,'.png'])
        export_fig(h(1),'-transparent','-r300',['.\fig\',figName,'.png'])
        h(1).Renderer='painters';
        print(h(1),'-depsc',['.\fig\',figName,'.eps'])
        cell2Save={'PaperType','PaperUnits','PaperSize','PaperPosition'};
        
        
        %print(h(1),'-dpdf',['.\fig\',figName,'.pdf'])
        export_fig(h(1),'-transparent',['.\fig\',figName,'.pdf'])
        
        h(1).Color=[0.93 0.93 0.93];
        for jj=1:numel(axh)
            axh(jj).Color=hCol{jj};
        end
        hgsave(h(1),['.\fig\',figName,'.fig'])
    end
    
    
end
