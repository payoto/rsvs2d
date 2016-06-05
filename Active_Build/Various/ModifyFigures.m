






%% Modify figures for postreatment
function []=ModifyFigures(h)

    ax=findobj(h,'type','axes');
    axes(ax(1));
    box=axis;
    
    for ii=1:length(ax)
       
        ax(ii).Visible='off';
        axes(ax(ii));
        axis equal
        axis(box)
    end
    
    
end


