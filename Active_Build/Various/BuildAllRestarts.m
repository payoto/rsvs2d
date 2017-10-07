function []=BuildAllRestarts(rootDir)
    
    [pathDir,nameDir]=FindDir(rootDir,'Dir',1);
    tot=0;
    totAll=0;
    for ii=1:numel(pathDir)
        
        [kk,kkf]=ChangeDirOptimParam(pathDir{ii},rootDir);
        tot=tot+kk;
        totAll=totAll+kkf+kk;
        disp(pathDir{ii})
        fprintf('%i of %i succesful\n',kk,kkf+kk)
        
    end

    fprintf('FINAL %i of %i succesful\n',tot,totAll)
    
    
end

function [kk,kkf]=ChangeDirOptimParam(pathDir,rootDir)
    
    [pathOptim,nameOptim]=FindDir(pathDir,'OptimRes',0);
    [pathParam]=regexprep(pathOptim,'OptimRes','FinalParam');
    restartPath=[pathDir,filesep,'iteration_0',filesep,'profile_0',filesep,'restart_0_0.mat'];
    kk=0;
    kkf=0;
    for ii=1:numel(pathOptim)
        try
            optimstruct=load(pathOptim{ii},'optimstruct');
            [optimstruct]=RewriteOptimForPC(optimstruct.optimstruct,rootDir);
            paramoptim=load(pathParam{ii},'paramoptim');
            [paramoptim]=RewriteParamForPC(paramoptim.paramoptim,rootDir);
            grid=load(restartPath,'grid');
            grid=grid.grid;
            newResPath=regexprep(pathOptim{ii},'OptimRes','RestartOptim');
            supportOptim=paramoptim.optim.supportOptim;
            [outinfo]=ReconstructOutinfo(optimstruct);
            save(newResPath,'optimstruct','paramoptim','supportOptim','grid','outinfo')
            kk=kk+1;
        catch MEid
            MEid.getReport
            kkf=kkf+1;
        end
    end
    
    
end

function [paramoptim]=RewriteParamForPC(paramoptim,rootDir)
    
    
    newCFDFolder=regexprep(ExtractVariables({'CFDfolder'},paramoptim),...
        '^.*source[_,a-Z]*',[regexprep(cd,'\\','\\\\'),filesep]);
    newresultRoot=regexprep(ExtractVariables({'resultRoot'},paramoptim.parametrisation),...
        '^.*source[_,a-Z]*',[regexprep(cd,'\\','\\\\'),filesep]);
    paramoptim=SetVariables({'CFDfolder'},{newCFDFolder},paramoptim);
    paramoptim.parametrisation=SetVariables({'resultRoot'},{newresultRoot},paramoptim.parametrisation);
    %'resultRoot', paramoptim.parametrisation
end