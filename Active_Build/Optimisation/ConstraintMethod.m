%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%            Optimisation Using
%          Parametric Snakes for
%         for Aerodynamic shape
%         parametrisation
%             Constraint Handling
%
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [population]=ConstraintMethod(entryPoint,paramoptim,population,varargin)
    % Function distributing the optimisation to various methods
    
    
    switch entryPoint
        case 'DesVar'
            varExtract={'desVarConstr','desVarVal'};
            [desVarConstr,desVarVal]=ExtractVariables(varExtract,paramoptim);
            
            for ii=1:length(desVarConstr)
                
                [population]=DesignVariableConsCaller(desVarConstr{ii},desVarVal{ii},...
                    paramoptim,population,varargin{:});
            end
            
        case 'Res'
            varExtract={'resConstr','resVal'};
            [resConstr,resVal]=ExtractVariables(varExtract,paramoptim);
            
            for ii=1:length(resConstr)
                
            end
    end
    
end

%% Design Variable cases

function [population]=DesignVariableConsCaller(constrName,constrVal,paroptim,population,varargin)
    
    switch constrName
        case 'MeanVolFrac'
            [population]=MeanVolumeFraction(constrVal,paroptim,population);
        case 'SumVolFrac'
            
        case 'Naca0012'
            
        otherwise
            error('Design Variable Constraint Not Recognised')
    end
            
    
    
end


function [population]=MeanVolumeFraction(constrVal,paroptim,population)
    
    varExtract={'desVarRange'};
    [desVarRange]=ExtractVariables(varExtract,paroptim);
    
    
    for ii=1:length(population)
       
        fillStart=population(ii).fill;
        
        meanFill=mean(fillStart);
        ratio=constrVal/meanFill;
        if ratio<=1
            population(ii).fill=fillStart*ratio;
        else
            maxFill=max(fillStart);
            if maxFill*ratio<=max(desVarRange)
                population(ii).fill=fillStart*ratio;
            else
                
                [population(ii).fill,population(ii).constraint]=...
                    IterativeMeanFill(fillStart,desVarRange,constrVal);
                
            end
        end
        
        
    end
    
    
    
end

function [fill,isConstr]=IterativeMeanFill(fill,desVarRange,constrVal)
    isConstr=true;
    ratio=mean(fill);
    kk=0;
    n=length(fill);
    while ratio~=constrVal && kk<=n+1;
        maxFill=max(desVarRange);
        
        fillBound=fill*ratio>=maxFill;
        fill(fillBound)=maxFill;
        fill(~fillBound)=fill(~fillBound)*ratio;
        ratio=mean(fill);
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
end



%% Results cases







