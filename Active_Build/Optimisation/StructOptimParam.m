%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%       Optimisation Using 
%          Parametric Snakes for
%      for Aerodynamic shape
%         parametrisation
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [paroptim]=StructOptimParam(caseStr)
    % Main function that allows changes
    
    [paroptim]=eval(caseStr);
    paroptim.general.case=caseStr;
    
    [paroptim.structdat]=ExploreStructureTree(paroptim);
    paroptim.structdat.vardat.names=[paroptim.structdat.vars(:).name];
    paroptim.structdat.vardat.varmatch=zeros(size(paroptim.structdat.vardat.names));
    for ii=1:length(paroptim.structdat.vars)
        jj=regexp(paroptim.structdat.vardat.names,paroptim.structdat.vars(ii).name);
        paroptim.structdat.vardat.varmatch(jj)=ii;
    end
    
    
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    paroptim.parametrisation=structInputVar(paramCase);
    
end


%% Deafult Optimisation Inputs

function [paroptim]=DefaultOptim()
    
    paroptim.general=DefaultOptimGeneral();
    
end

function [paroptimgeneral]=DefaultOptimGeneral()
    
    paroptimgeneral.paramCase='optimDefault';
    paroptimgeneral.optimMethod='DE';
    
end




%% Standard Modifications



%% Callable functions

function [paroptim]=StandardOptim()
    
    [paroptim]=DefaultOptim();
    
end

