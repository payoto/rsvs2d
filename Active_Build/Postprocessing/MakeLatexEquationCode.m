function [cellStr]=MakeLatexEquationCode(symIn,pre,post)
    
    cellStr=cell([1,numel(symIn)]);
    for ii=1:numel(symIn)
        str='';
        str=char(str,'\begin{equation}');
        if numel(pre)>0
            str=char(str,pre);
        end
        str=char(str,latex(simplify(symIn(ii))));
        if numel(post)>0
            str=char(str,post);
        end
        str=char(str,'\end{equation}');
        cellStr{ii}=str;
    end
    
    
end
