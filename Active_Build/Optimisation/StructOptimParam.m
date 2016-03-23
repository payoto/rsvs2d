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
    paroptim.general.optimCase=caseStr;

end

function structdat=GetStructureData(paroptim)
        
    [structdat]=ExploreStructureTree(paroptim);
    structdat.vardat.names=[structdat.vars(:).name];
    structdat.vardat.varmatch=zeros(size(structdat.vardat.names));
    for ii=1:length(structdat.vars)
        jj=regexp(structdat.vardat.names,structdat.vars(ii).name);
        structdat.vardat.varmatch(jj)=ii;
    end
    
end
%% Deafult Optimisation Inputs

function [paroptim]=DefaultOptim()
    
    paroptim.general=DefaultOptimGeneral();
    paroptim.structdat=GetStructureData(paroptim);
    
end

function [paroptimgeneral]=DefaultOptimGeneral()
    
    paroptimgeneral.optimCase='optimDefault';
    paroptimgeneral.paramCase='optimDefault';
    paroptimgeneral.optimMethod='DE';
    paroptimgeneral.desVarRange=[0,1];
    paroptimgeneral.nDesVar=[0];
    paroptimgeneral.nPop=[10];
    paroptimgeneral.maxIter=10;
    
end



%% Standard Modifications

function [paraminit]=ChangeSnakeInit(paraminit)
    
    paraminit.snakes.step.snakesSteps=20;
    paraminit.general.restart=false;
end


%% Callable functions

function [paroptim]=StandardOptim()
    
    [paroptim]=DefaultOptim();
    
    varExtract={'paramCase',};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    
    paroptim.parametrisation=structInputVar(paramCase);
    paroptim.initparam=ChangeSnakeInit(paroptim.parametrisation);
    
end

