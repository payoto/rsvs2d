
function []=FigureTextTools(figList,varargin)
    
    [interpreter,FontSizeDelta]=HandleVarargin(varargin);
    LatexFigure(figList,interpreter);
end


function []=LatexFigure(figList,interpreter)
    % turns all the text to the latex interpreter
    if nargin <2
        interpreter='latex';
    end
    
    for ii=figList;
        h=findobj(ii,'type','axes');
        if numel(h)>0
            for jj=1:numel(h)
                h(jj).TickLabelInterpreter=interpreter;
                h(jj).XLabel.Interpreter=interpreter;
                h(jj).YLabel.Interpreter=interpreter;
                try
                    h(jj).ZLabel.Interpreter=interpreter;
                catch
                end
            end
        end
        h=findobj(ii,'type','legend');
        if numel(h)>0
            for jj=1:numel(h)
                h(jj).Interpreter=interpreter;
            end
        end
        h=findobj(ii,'type','text');
        if numel(h)>0
            for jj=1:numel(h)
                h(jj).Interpreter=interpreter;
            end
        end
    end
end


function [interpreter,FontSizeDelta]=HandleVarargin(cellArgin)
    
    interpreter='latex';
    FontSizeDelta=0;
        
    for ii=1:2:numel(cellArgin)
        eval([cellArgin{ii},'=cellArgin{ii+1};']);
    end
    
end
