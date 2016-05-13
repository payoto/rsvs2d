function [] = include_Optimisation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%%

function [desVarList]=ExtractActiveVariable(nFill,notDesInd,inactiveVar)
    
    desVarList=1:nFill;
    desVarList([notDesInd,inactiveVar])=[];
    
end




function [inactiveVar]=SelectInactiveVariables(newFill,varActive)
    
    switch varActive
        case 'all'
            inactiveVar=[];
        case 'border'
            is0=newFill==0;
            is1=newFill==1;
            
            inactiveVar=find(is0);
            inactiveVar=[inactiveVar,find(is1)];
            
        case 'wideborder'
            error('Not coded yet')
            
        otherwise
            error('unrecognised variable activation criterion')
    end
    
    
end
