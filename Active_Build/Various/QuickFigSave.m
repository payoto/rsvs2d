function []=QuickFigSave(figList)
    
    for ii=figList;
        h=findobj(ii);
        h(1).Renderer='opengl';
        h(1).Color='none';
        hCol=h(3).Color;
        h(3).Color='none';
        figName=matlab.lang.makeValidName([h(1).Name,int2str(ii)]);
        print(h(1),'-r300','-dpng',['.\fig\',figName,'.png'])
        h(1).Renderer='painters';
        print(h(1),'-depsc',['.\fig\',figName,'.eps'])
        h(1).Color=[0.93 0.93 0.93];
        h(3).Color=hCol;
        hgsave(h(1),['.\fig\',figName,'.fig'])
    end
    
    
end