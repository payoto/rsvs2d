function []=ModifyAllOptimStruct(rootDir)
    
    [pathDir,nameDir]=FindDir(rootDir,'Dir',1);
    
    for ii=1:numel(pathDir)
        
        ChangeDirOptimParam(pathDir{ii},rootDir);
    end
    
    
    
end

function []=ChangeDirOptimParam(pathDir,rootDir)
    
    [pathOptim,nameOptim]=FindDir(pathDir,'OptimRes',0);
    
    for ii=1:numel(pathOptim)
        instruct=load(pathOptim{ii});
        [optimstruct]=RewriteOptimForPC(instruct.optimstruct,rootDir);
        save(pathOptim{ii},'optimstruct')
        
    end
    
    [pathParam,nameParam]=FindDir(pathDir,'FinalParam',0);
    
    for ii=1:numel(pathOptim)
        instruct=load(pathParam{ii});
        [paramoptim]=RewriteParamForPC(instruct.paramoptim,rootDir);
        save(pathParam{ii},'paramoptim')
        
    end
    
end

function [paramoptim]=RewriteParamForPC(paramoptim,rootDir)
    
    
    newCFDFolder=regexprep(ExtractVariables({'CFDfolder'},paramoptim),...
        '^.*source',cd);
    paramoptim=SetVariables({'CFDfolder'},{newCFDFolder},paramoptim);
    %'resultRoot', paramoptim.parametrisation
end