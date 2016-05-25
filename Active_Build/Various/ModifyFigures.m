






%% Modify figures for postreatment
function []=ModifyFigures(h)

    ax=findobj(h,'type','axes');
    
    for ii=1:length(ax)
       
        set(ax(ii),'fontsize',12)
        set(ax(ii),'fontsize',12)
        ax(ii).YLabel.FontSize = 14;
        ax(ii).XLabel.FontSize = 14;
        ax(ii).ZLabel.FontSize = 14;
    end
    
    
end


