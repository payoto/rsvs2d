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
    
    % Error Report
    [FIDError]=OpenErrorReportFile(writeDirectory,marker);
    GenerateErrorReportFile(t,marker,FIDError)
end

function [writeDirectory]=OptimisationOutput_profile(out,nIter,nProf,loop,restartsnak,snakSave)
    
    marker=out.marker;
    t=out.tOutput;
    rootDir=out.rootDir;
%     iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yymmddTHHMM')];
%     profStr=['\profile_',int2str(nProf),'_',datestr(t,'yymmddTHHMMSS')];
    markerShort=[int2str(nIter),'_',int2str(nProf)];
    
    iterStr=['\iteration_',int2str(nIter)];
    profStr=['\profile_',int2str(nProf)];
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
    WriteFullInfoProfile(writeDirectory,nProf,marker,t,nIter)
end

function [out]=OptimisationOutput_iteration(nIter,out,population,errorReports)
    
    t=out.tOutput;
    rootDir=out.rootDir;
%     marker=['iteration_',int2str(nIter),datestr(t,'_yymmddTHHMM')];
%     iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yymmddTHHMM')];
    marker=['iteration_',int2str(nIter)];
    iterStr=['\iteration_',int2str(nIter)];
    writeDirectory=[rootDir,iterStr];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    % iter entry
    [FIDIter]=OpenIterIndexFile(rootDir,out.marker);
    returnEntries=GenerateIterIndexEntry(FIDIter,nIter,population);
    % Error Info
    [FIDError]=OpenErrorReportFile(rootDir,out.marker);
    GenerateErrorReportEntries(FIDError,nIter,errorReports,returnEntries)
    
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
    
    paramCell{1}='# Iteration Index File';
    paramCell{2}=['# ',datestr(t)];
    paramCell{3}=['# ',marker];
    paramCell{4}=[' '];
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function [returnEntries]=GenerateIterIndexEntry(FID,nIter,population)
    
    fieldsAdd=fieldnames(population(1).additional);
    nAdditional=length(fieldsAdd);
    entryCell{1+length(population)}={};
    entryCell{1}=['# ',datestr(now),', iter , member , objective , constraint '];
    for ii=1:nAdditional
        entryCell{1}=[entryCell{1},' , ',fieldsAdd{ii}];
    end
    entryCell{1}=[entryCell{1},', exception , fill '];
    
    for ii=1:length(population)
        
        str=[int2str(nIter)];
        str=[str,' , ', int2str(ii)];
        str=[str,' , ', num2str(population(ii).objective,' %12.7e ')];
        str=[str,' , ', int2str(population(ii).constraint)];
        
        for jj=1:nAdditional
            str=[str,' , ', num2str(population(ii).additional.(fieldsAdd{jj}),' %12.7e ')];
        end
        str=[str,' , ', population(ii).exception];
        str=[str,' , ', num2str(population(ii).fill,' %12.7e ')];
        entryCell{ii+1}=str;
    end
    
    WriteToFile(entryCell,FID);
    fclose(FID);
    returnEntries=entryCell(2:end);
end

function [FID]=OpenIterIndexFile(rootOptim,marker)
    % Creates a file in the current directory to write data to.
    
    rootOptim=MakePathCompliant(rootOptim);
    fileName=['Index_',marker,'.txt'];
    FID=fopen([rootOptim,filesep,fileName],'a');
    
end

function [FID]=OpenErrorReportFile(rootOptim,marker)
    % Creates a file in the current directory to write data to.
    
    rootOptim=MakePathCompliant(rootOptim);
    fileName=['ErrorReport_',marker,'.txt'];
    FID=fopen([rootOptim,filesep,fileName],'a');
    
end

function []=GenerateErrorReportFile(t,marker,FID)
    
    paramCell{1}='# Error Report File';
    paramCell{2}=['# ',datestr(t)];
    paramCell{3}=['# ',marker];
    paramCell{4}=[' '];
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function []=GenerateErrorReportEntries(fID,nIter,errorReports,indexEntries)
    
    writeReport{1}='------------------------------------------------------------';
    writeReport{2}=['   ITERATION ',int2str(nIter)];
    writeReport{3}='------------------------------------------------------------';
    kk=4;
    for ii=1:length(errorReports)
        
        if ~isempty(errorReports{ii})
            writeReport{kk}=indexEntries{ii};
            kk=kk+1;
            writeReport{kk}=errorReports{ii};
            kk=kk+1;
        end
        
    end

    WriteToFile(writeReport,fID);
    fclose(fID);
    
end

function []=WriteFullInfoProfile(writeDirectory,nProf,marker,t,nIter)
    fid=fopen([writeDirectory,filesep,'InfoProfile_',int2str(nIter),'_',int2str(nProf),'.dat'],'w');
    
    ii=1;
    cellDat{ii}=marker;ii=ii+1;
    cellDat{ii}=datestr(t);ii=ii+1;
    cellDat{ii}=['iteration : ',int2str(nIter)]; ii=ii+1;
    cellDat{ii}=['Profile : ',int2str(nProf)]; ii=ii+1;
    
    WriteToFile(cellDat,fid);
    fclose(fid);
end

%% 
function [h]=OptimHistory(optimstruct,knownOptim,dirOptim)
    
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1)
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar]);
    hold on
    for ii=1:nIter
        iterRes(ii,:)=[optimstruct(ii).population(:).objective];
        lSub1(1)=plot(ones(1,nVar)*ii,iterRes(ii,:),'b*');
    end
    switch dirOptim
        case 'min'
            minRes=min(iterRes,[],2);
        case 'max'
            minRes=max(iterRes,[],2);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:nIter,minRes,'r-');
    lSub1(3)=plot(1:nIter,meanRes,'color',[0.7 0 0]);
    
    
    lSub1(4)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Population ',dirOptim,'imum'],'Population mean',...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex')
    ylabel('$J(\mathbf{x})$', 'interpreter','latex')
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,nIter+1];
    
    
    axis(box);
    
    % Plot 2
    subplot(1,2,2)
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            errMean=-(knownOptim-meanRes);
        case 'max'
            errMeasure=knownOptim-minRes;
            errMean=(knownOptim-meanRes);
    end
    
    lSub2(1)=semilogy(1:nIter,errMeasure);
    hold on
    lSub2(2)=semilogy(1:nIter,errMean);
    
    for ii=1:nIter
        stdline=[errMean(ii)-stdRes(ii),errMean(ii)+stdRes(ii)];
        lSub2(3)=semilogy([ii ii],stdline,'g--+');
    end
    legend(lSub2,{['Population ',dirOptim,'imum'],'Population mean',...
        'Standard Deviation'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex')
    ylabel('$J^*_T-J(\mathbf{x})$', 'interpreter','latex')
    
    
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
