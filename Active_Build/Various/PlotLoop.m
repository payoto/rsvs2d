function [varargout]=PlotLoop(loop,typeLoop,isfill)
    
    if nargin<3
       isfill=false; 
    end
    
    colors = get(gca,'ColorOrder');
    index  = get(gca,'ColorOrderIndex');
    hold on
    n_colors = size(colors,1);
    if index > n_colors
        index = 1;
    end
    next_color = colors(index,:);
    if ~isfill
        for ii=1:numel(loop)
            l(ii)=plot(loop(ii).(typeLoop)([1:end,1],1),loop(ii).(typeLoop)([1:end,1],2),'Color',next_color);
        end
    else
        for ii=1:numel(loop)
            l(ii)=patch(loop(ii).(typeLoop)([1:end,1],1),loop(ii).(typeLoop)([1:end,1],2),next_color);
        end
    end
    set(gca,'ColorOrderIndex',index+1);
    if nargout>0
        varargout{1}=l;
    end
end