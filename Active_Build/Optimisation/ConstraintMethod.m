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
                [population]=ResultVariableConsCaller(resConstr{ii},resVal{ii},...
                    paramoptim,population,varargin{:});
            end
    end
    
end

%% Design Variable cases

function [population]=DesignVariableConsCaller(constrName,constrVal,paroptim,population,varargin)
    
    switch constrName
        case 'MeanVolFrac'
            [population]=MeanVolumeFraction(constrVal,paroptim,population);
        case 'MinSumVolFrac'
            [population]=MinSumVolumeFraction(constrVal,paroptim,population);
        case 'Naca0012'
        
        case ' '
            
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
    ratio=constrVal/mean(fill);
    kk=0;
    n=length(fill);
    while ratio~=constrVal && kk<=n+1;
        maxFill=max(desVarRange);
        
        fillBound=((fill*ratio)>=maxFill);
        fill(fillBound)=maxFill;
        fill(~fillBound)=fill(~fillBound)*ratio;
        ratio=constrVal/mean(fill);
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
end


function [population]=MinSumVolumeFraction(constrVal,paroptim,population)
    
    varExtract={'desVarRange'};
    [desVarRange]=ExtractVariables(varExtract,paroptim);
    
    
    for ii=1:length(population)
       
        fillStart=population(ii).fill;
        
        sumFill=sum(fillStart);
        ratio=constrVal/sumFill;
        if ~isfinite(ratio)
            population(ii).fill=fillStart;
            population(ii).constraint=false;
        elseif ratio<=1
            population(ii).fill=fillStart;
        elseif ratio>1
            maxFill=max(fillStart);
            if maxFill*ratio<=max(desVarRange)
                population(ii).fill=fillStart*ratio;
            else
                
                [population(ii).fill,population(ii).constraint]=...
                    IterativeMinFill(fillStart,desVarRange,constrVal);
                
            end
    
        end
        
        
    end
    
    
    
end
function [fill,isConstr]=IterativeMinFill(fill,desVarRange,constrVal)
    isConstr=true;
    ratio=constrVal/sum(fill);
    kk=0;
    n=length(fill);
    while ratio>1 && kk<=n+1;
        maxFill=max(desVarRange);
        
        fillBound=((fill*ratio)>=maxFill);
        fill(fillBound)=maxFill;
        fill(~fillBound)=fill(~fillBound)*ratio;
        ratio=constrVal/sum(fill);
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
end
%% Results cases

function [population]=ResultVariableConsCaller(constrName,constrVal,paroptim,population,varargin)
    
    switch constrName
        case 'AeroResidual'
            population=CheckAerodynamicResidual(constrVal,population);
        case 'AeroResidualBarrier'
            population=BarrierAerodynamicResidual(constrVal,population);
        case 'AeroLift'
            
        case 'Volume'
        
        case ' '
            
        otherwise
            error('Design Variable Constraint Not Recognised')
    end
            
    
    
end

function population=CheckAerodynamicResidual(constrVal,population)
    
    for ii=1:length(population)
        
        if population(ii).constraint
           constrViolation= (population(ii).additional.res>constrVal) && ...
               (population(ii).additional.res~=0);
           population(ii).constraint=~constrViolation;
        end
        
    end
    
    
end

function population=BarrierAerodynamicResidual(constrVal,population)
    
    for ii=1:length(population)
        
        if population(ii).constraint>0
            if population(ii).additional.res<constrVal(1) ... % Constraint fully satisfied
                    || population(ii).additional.res==0
                population(ii).constraint=1;

            elseif population(ii).additional.res>constrVal(2) % Constraint fully violated
                population(ii).constraint=0;

            else
                population(ii).constraint=1-(0.5+0.5*tanh(3/(-(constrVal(1)-constrVal(2))/2)...
                    *population(ii).additional.res));
            end
        end
        
    end
    
    
end



