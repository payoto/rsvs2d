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
        case 'final'
            out=OptimisationOutput_Final(paramoptim,varargin{:});
    end
    
    
end


%% 

function [marker,t, writeDirectory]=OptimisationOutput_Init(paramoptim)
    
    % Unpack necessary variables
    varExtract={'optimCase','typDat','resultRoot','archiveName'};
    [optimCase]=ExtractVariables(varExtract(1),paramoptim);
    varExtract={'optimCase','typDat','resultRoot','archiveName'};
    [typDat,resultRoot,archiveName]=ExtractVariables(varExtract(2:end),paramoptim.parametrisation);
    
    [marker,t]=GenerateResultMarker([optimCase]);
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
    
    % Iter Index
    [FIDIter]=OpenIterIndexFile(writeDirectory,marker);
    GenerateIterIndexFile(t,marker,FIDIter)
end

function [writeDirectory]=OptimisationOutput_profile(out,nIter,nProf,loop,restartsnak,snakSave)
    
    marker=out.marker;
    t=out.tOutput;
    rootDir=out.rootDir;
    iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yymmddTHHMM')];
    profStr=['\profile_',int2str(nProf),'_',datestr(t,'yymmddTHHMMSS')];
    markerShort=[int2str(nIter),'_',int2str(nProf)];
    writeDirectory=[rootDir,iterStr,profStr];
    writeDirectory=MakePathCompliant(writeDirectory);
    mkdir(writeDirectory);
    
    savStruct.restartsnak=restartsnak;
    savStruct.snakSave=snakSave;
    savStruct.loop=loop;
    
     % Output boundary data file
    [fidBoundary]=OpenBoundaryFile(writeDirectory,markerShort);
    for ii=1:length(loop)
        loop(ii).subdivision=loop(ii).subdivspline;
    end
    BoundaryOutput(loop,fidBoundary);
    fclose(fidBoundary);
    
    
    GenerateProfileBinary(writeDirectory,markerShort,savStruct)
end

function [out]=OptimisationOutput_iteration(nIter,out,population)
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=['iteration_',int2str(nIter),datestr(t,'_yymmddTHHMM')];
    iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yymmddTHHMM')];
    writeDirectory=[rootDir,iterStr];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    % iter entry
    [FIDIter]=OpenIterIndexFile(rootDir,out.marker);
    GenerateIterIndexEntry(FIDIter,nIter,population);
    
    CopyDiary(writeDirectory,marker)
    GeneratePopulationBinary(writeDirectory,marker,population)
    h=CheckOptimProfile('iter_all',writeDirectory);
   %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
   figName=[writeDirectory,'\profiles_',marker,'.fig'];
    figName=MakePathCompliant(figName);
   hgsave(h,figName);
    close(h);
    
    
end

function [out]=OptimisationOutput_Final(paroptim,out,optimstruct)
    
    varExtract={'direction','knownOptim'};
    [direction,knownOptim]=ExtractVariables(varExtract,paroptim);
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=out.marker;
    markerSmall=datestr(t,'_yymmddTHHMM');
    writeDirectory=[rootDir];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    CopyDiary(writeDirectory,marker)
    GenerateIterResultBinary(writeDirectory,marker,optimstruct)
    GenerateOptimalSolDir(writeDirectory,markerSmall,direction,optimstruct)
    [h]=OptimHistory(optimstruct,knownOptim,direction);
    %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
    figName=[writeDirectory,'\Optimisation_',marker,'.fig'];
    figName=MakePathCompliant(figName);
    hgsave(h,figName);
   
end


%% Index File ops


function []=GenerateIterIndexFile(t,marker,FID)
    
    paramCell{1}='# Parameter File';
    paramCell{2}=['# ',datestr(t)];
    paramCell{3}=['# ',marker];
    paramCell{4}=[' '];
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function []=GenerateIterIndexEntry(FID,nIter,population)
    
    fieldsAdd=fieldnames(population(1).additional);
    nAdditional=length(fieldsAdd);
    entryCell{length(population)}={};
    entryCell{1}=['# ',datestr(now),', iter , member , objective , constraint, fill'];
    for ii=1:nAdditional
        entryCell{1}=[entryCell{1},' , ',fieldsAdd{ii}];
    end
    
    for ii=1:length(population)
        
        str=[int2str(nIter)];
        str=[str,' , ', int2str(ii)];
        str=[str,' , ', num2str(population(ii).objective,' %12.7e ')];
        str=[str,' , ', int2str(population(ii).constraint)];
        str=[str,' , ', num2str(population(ii).fill,' %12.7e ')];
        
        for jj=1:nAdditional
            str=[str,' , ', num2str(population(ii).additional.(fieldsAdd{jj}),' %12.7e ')];
        end
        
        entryCell{ii+1}=str;
    end
    
    WriteToFile(entryCell,FID);
    fclose(FID);
end

function [FID]=OpenIterIndexFile(rootOptim,marker)
    % Creates a file in the current directory to write data to.
    
    rootOptim=MakePathCompliant(rootOptim);
    fileName=['Index_',marker,'.txt'];
    FID=fopen([rootOptim,filesep,fileName],'a');
    
end

%% 

function [h]=OptimHistory(optimstruct,knownOptim,dirOptim)
    
    
    h=figure('Name','Optimisation Result','Position',[20 100 1300 800]);
    subplot(1,2,1)
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar]);
    hold on
    for ii=1:nIter
        iterRes(ii,:)=[optimstruct(ii).population(:).objective];
        plot(ones(1,nVar)*ii,iterRes(ii,:),'b*')
    end
    minRes=min(iterRes,[],2);
    plot(1:nIter,minRes,'r-')
    
    
    plot([0,nIter],[knownOptim knownOptim],'r--')
    subplot(1,2,2)
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
        case 'max'
            errMeasure=knownOptim-minRes;
    end
    semilogy(1:nIter,errMeasure)
    
end

function []=GenerateProfileBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,'\restart_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'-struct','restartstruct');
    
end

function []=GeneratePopulationBinary(resultDirectory,marker,population)
    
    fileName=[resultDirectory,'\population_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'population');
    
end

function []=GenerateIterResultBinary(resultDirectory,marker,optimstruct)
    
    fileName=[resultDirectory,'\OptimRes_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'optimstruct');
    
end

function []=GenerateOptimalSolDir(resultDirectory,markerSmall,optimDirection,optimstruct)
    
    [~,posOpt]=eval([optimDirection,'([optimstruct(end).population(:).objective])']);
    optimsolution=optimstruct(end).population(posOpt);
    resultDirectory=[resultDirectory,'\Optimal_',markerSmall];
    resultDirectory=MakePathCompliant(resultDirectory);
    profileDir=optimsolution.location;
    
    copyfile(profileDir,resultDirectory);
    
    c=dir(resultDirectory);
    isFileName=false;
    ii=0;
    while(~isFileName)
        ii=ii+1;
        isFileName=~isempty(regexp(c(ii).name,'restart', 'once'));
    end
    load([resultDirectory,filesep,c(ii).name])
    h=CheckOptimProfile('loop',loop);
    hgsave(h,MakePathCompliant([resultDirectory,'\OptProf_',markerSmall,'.fig']));
    fileName=[resultDirectory,'\OptProf_',markerSmall,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'optimsolution');
    
end
