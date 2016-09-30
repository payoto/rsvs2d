%% Refined Optimisation steps

function []=HandleRefinement(paramoptim,iterstruct)
    
    % grid refinement
    
    % parameter update
    
    % grid matching
    
    % Fill matching
    
    
end


% 


function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        InitialiseRefinement(paramoptim,oldBase,refStep)
    
    procStr='REFINE OPTIMISATION PROCESS';
    [tStart]=PrintStart(procStr,1);
    
    % Initialise Optimisation
    % Get Parametrisation parameters
%     [outinfo]=OptimisationOutput('init',paramoptim);
    warning ('[~,paramoptim]=ConstraintMethod(''init'',paramoptim,[]); Not supported');

    % Initialise Grid
    varNames={'refineCellLvl'};
    refineCellLvl=ExtractVariables(varNames,paroptim.parameterisation);
    refparamsnake=SetVariables({'refineGrid'},{refineCellLvl},paroptim.parameterisation);
    
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(refparamsnake);
    
    newgrid.base=baseGrid;
    newgrid.refined=gridrefined;
    newgrid.connec=connectstructinfo;
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(refparamsnake,gridrefined);
    
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [paramoptim.general.desvarconnec]=...
        ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    
    
    
    
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim);
    
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,...
        paramoptim.spline,outinfo,0,0,0);
    
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    
    [~]=PrintEnd(procStr,1,tStart);
end