function []=CompareParams(param1,param2)
    try
        vars={param1.structdat.vars.name};
        isparam=true;
    catch
        vars=fieldnames(param1);
        isparam=false;
    end
    isEqualStruct=false(size(vars));
    for ii=1:numel(vars)
        if isparam
            comp1=ExtractVariables(vars(ii),param1);
            comp2=ExtractVariables(vars(ii),param2);
        else
            comp1=param1.(vars{ii});
            comp2=param2.(vars{ii});
            
        end
        
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
        %isequal=true;
        isequal=numel(comp1)==numel(comp2);
        if isequal
            for ii=1:numel(comp1)
                [isequal(ii)]=TestEqual(comp1{ii},comp2{ii});
            end
        end
        isequal=all(isequal);
    elseif isstruct(comp1)
        vars=fieldnames(comp1);
        isequalarray=numel(comp1)==numel(comp2);
        if isequalarray
            for ii=1:numel(comp1)
                for jj=1:numel(vars)
                    [isequalarray(ii,jj)]=TestEqual(comp1(ii).(vars{jj}),comp2(ii).(vars{jj}));
                end
                
            end
        end
        
        isequal=all(all(isequalarray));
    end
    
    
end