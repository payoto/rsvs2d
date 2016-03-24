%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          snake parametrisation
%      for Aerodynamic shape parametrisation
%           - Outputs Management Function -
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [out]=OptimisationOutput(entryPoint,paramoptim,varargin)
    
    
    switch entryPoint
        case 'init'
            [out.marker,out.tOutput,out.rootDir]=OptimisationOutput_Init(paramoptim);
        case 'profile'
            out=varargin{1};
            out.dirprofile=OptimisationOutput_profile(varargin{:});
        case 'iteration'
            out=OptimisationOutput_iteration(varargin{:});
    end
    
    
end


%% 

function [marker,t, writeDirectory]=OptimisationOutput_Init(paramoptim)
    
    % Unpack necessary variables
    varExtract={'optimCase','typDat','resultRoot','archiveName'};
    [optimCase]=ExtractVariables(varExtract(1),paramoptim);
    varExtract={'optimCase','typDat','resultRoot','archiveName'};
    [typDat,resultRoot,archiveName]=ExtractVariables(varExtract(2:end),paramoptim.parametrisation);
    
    [marker,t]=GenerateResultMarker([optimCase,'_',typDat]);
    [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t);
    
    % Parameter Data
    [fidParam]=OpenParamFile(writeDirectory,marker);
    GenerateParameterFile(fidParam,paramoptim,t,marker);
    
    [fidParam]=OpenParamFile(writeDirectory,[marker,'_parametrisation']);
    GenerateParameterFile(fidParam,paramoptim.parametrisation,t,[marker,'_parametrisation']);
    
    % Comments file
    [fidComments]=OpenCommentsFile(writeDirectory,marker);
    indexEntry=MakeCommentsFile(fidComments,paramoptim.parametrisation,t,writeDirectory);
    fclose(fidComments);
    
    % Index File Entry 
    [fidIndexFile]=OpenIndexFile(resultRoot,archiveName);
    WriteToFile(indexEntry,fidIndexFile);
    fclose(fidIndexFile);
end

function [writeDirectory]=OptimisationOutput_profile(out,nIter,nProf,loop,restartsnak,snakSave)
    
    marker=out.marker;
    t=out.tOutput;
    rootDir=out.rootDir;
    iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yyyy-mm-ddTHHMM')];
    profStr=['\profile_',int2str(nProf),'_',datestr(t,'yyyy-mm-ddTHHMMSS')];
    writeDirectory=[rootDir,iterStr,profStr];
    system(['md "',writeDirectory,'"']);
    
     % Output boundary data file
    [fidBoundary]=OpenBoundaryFile(writeDirectory,marker);
    BoundaryOutput(loop,fidBoundary);
    fclose(fidBoundary);
    
    savStruct.restartsnak=restartsnak;
    savStruct.snakSave=snakSave;
    savStruct.loop=loop;
    GenerateProfileBinary(writeDirectory,marker,savStruct)
end


function [out]=OptimisationOutput_iteration(nIter,out,population)
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=['iteration_',int2str(nIter),datestr(t,'_yyyy-mm-ddTHHMM')];
    iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yyyy-mm-ddTHHMM')];
    writeDirectory=[rootDir,iterStr];
    
    CopyDiary(writeDirectory,marker)
    GeneratePopulationBinary(writeDirectory,marker,population)
    h=CheckOptimProfile('iter_all',writeDirectory);
   %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
   hgsave(h,[writeDirectory,'\profiles_',marker,'.fig']);
    close(h);
end

%% 

function []=GenerateProfileBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,'\restart_',marker,'.mat'];
    save(fileName,'-struct','restartstruct');
    
end

function []=GeneratePopulationBinary(resultDirectory,marker,population)
    
    fileName=[resultDirectory,'\population_',marker,'.mat'];
    save(fileName,'population');
    
end