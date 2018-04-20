function [ASOstruct]=ASOPerformanceAPI(optimstruct)
    
end

function [h]=PlotASOPerformance(ASOstruct,axDeOpt)
    
    h=figure('Name','ASO Performance');
    
    fieldsASO=fieldnames(ASOstruct);
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec','objFuncCalls','CD0'};
    %              scalar    vec       vec    vec     vec(numLvl)  vec(numLvl)   scalar     scalar  
    
    % Iteration plots overlaid on DEresults
    for ii=1:numel(ASOstruct)
        plot(axDeOpt,[ASOstruct(ii).DEIter,ASOstruct(ii).DEIter+1],...
            ASOstruct(ii).obj([1,end]),'k--');
        % Possibly change this to have a changing color with better
        % improvements
    end
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    % 2) Change in objective thanks to ASO (cloud)
    
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    % 4) Change of objective other the first iteration (normalised by total Delta)
    
    
    % #iterations needed to achieve 95% of the improvement vs Number of Desvar
    
    
    
    % Magnitude of the geometric step vs number of surf points and # of design variables
    
end