% Function for the postreatment of incomplete optimisation run
% The goal is to reconstruct enough information to take the run through 

function [outinfo,paramoptim,iterstruct]=PostTreatIncomplete(pathStr,nIter,iterstruct)
    include_Utilities
    include_PostProcessing
    include_Optimisation
    
    % Reconstruct outinfo
%     outinfo.rootDir=pathStr;
%     [outinfo.tOutput,outinfo.marker]=FindTime(pathStr);
    [outinfo]=ReconstructOutinfo(iterstruct);
    % iterstruct paramoptim
    [paramoptim2]=ReconstructParameter(pathStr,outinfo(1).marker);
    try 
        paramoptim=StructOptimParam(ExtractVariables({'optimCase'},paramoptim2));
        outVars={paramoptim2.structdat(:).vars.name};
        for ii=numel(outVars):-1:1
            [valVars{ii}]=ExtractVariables(outVars(ii),paramoptim2);
        end
        for ii=1:numel(outVars)
            try
            paramoptim=SetVariables(outVars(ii),valVars(ii),paramoptim);
            catch
            end
        end
    catch MEId
        disp(MEId.getReport)
        [paramoptim]=ReconstructParameter(pathStr,outinfo(1).marker);
    end
    % Reconstruct iterstruct
    if nargin<3
        [iterstruct]=ReconstructIterationStructure(pathStr,nIter,paramoptim);
    end
    
    
    OptimisationOutput('final',paramoptim,outinfo,iterstruct);
end


%% Standard Postreatment functions

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        warning('Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
      
    
    
end
