function []=CompareParams(param1,param2)
    
    vars={param1.structdat.vars.name};
    isEqualStruct=false(size(vars));
    for ii=1:numel(vars)
        comp1=ExtractVariables(vars(ii),param1);
        comp2=ExtractVariables(vars(ii),param2);
        
        isEqualStruct(ii)=TestEqual(comp1,comp2);
        
    end
    
    char(vars(~isEqualStruct))
    
end


function [isequal]=TestEqual(comp1,comp2)
    if isnumeric(comp1) || islogical(comp1)
        isequal=all(all(comp1==comp2));
    elseif ischar(comp1)
        isequal=strcmp(comp1,comp2);
    elseif iscell(comp1)
        isequal=true;
        for ii=1:numel(comp1)
            [isequal(ii)]=TestEqual(comp1{ii},comp2{ii});
        end
        isequal=all(isequal);
    else
        isequal=true;
    end
    
    
end