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
        case 'finalpost'
            out=OptimisationOutput_Final_Post(paramoptim,varargin{:});
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
    
    % Lay File
    [FID]=OpenOptimLayFile(writeDirectory,marker);
    pltFile=[writeDirectory,filesep,'Tec360PLT_',marker,'.plt'];
    PersnaliseLayFile(FID,pltFile);
end

function [writeDirectory]=OptimisationOutput_profile(out,nIter,nProf,loop,...
        restartsnak,snakSave,tecStruct)
    
    
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
    if nIter>0
        TecplotPortion_Profile(nIter,tecStruct.nPop,nProf,writeDirectory,tecStruct.baseGrid,...
            snakSave(end).volumefraction,restartsnak.snaxel,tecStruct.snakposition);
    else
        TecplotPortion_Init(nIter,tecStruct.nPop,nProf,writeDirectory,tecStruct.baseGrid,...
            tecStruct.fineGrid,snakSave(end).volumefraction,restartsnak.snaxel,...
            tecStruct.snakposition);
    end
    
    GenerateProfileBinary(writeDirectory,markerShort,savStruct)
    WriteFullInfoProfile(writeDirectory,nProf,marker,t,nIter)
end

function [out]=OptimisationOutput_iteration(nIter,out,population,errorReports)
    
    t=out.tOutput;
    rootDir=out.rootDir;
    fullMark=out.marker;
%     marker=['iteration_',int2str(nIter),datestr(t,'_yymmddTHHMM')];
%     iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yymmddTHHMM')];
    marker=['iteration_',int2str(nIter)];
    iterStr=['\iteration_',int2str(nIter)];
    writeDirectory=[rootDir,iterStr];
    writeDirectory=MakePathCompliant(writeDirectory);
    tecFilePath=[rootDir,filesep,'Tec360PLT_',fullMark,'.plt'];
    tecFilePath=MakePathCompliant(tecFilePath);
        ConcatenateTecplotFile(writeDirectory,tecFilePath)
        
    if nIter~=0
        % iter entry
        [FIDIter]=OpenIterIndexFile(rootDir,out.marker);
        returnEntries=GenerateIterIndexEntry(FIDIter,nIter,population);
        % Error Info
        [FIDError]=OpenErrorReportFile(rootDir,out.marker);
        GenerateErrorReportEntries(FIDError,nIter,errorReports,returnEntries)
        if ~isdir(writeDirectory)
            mkdir(writeDirectory);
        end
        CopyDiary(writeDirectory,marker)
        GeneratePopulationBinary(writeDirectory,marker,population)
        
        
        
        h=CheckOptimProfile('iter_all',writeDirectory);
        %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
        figName=[writeDirectory,'\profiles_',marker,'.fig'];
        figName=MakePathCompliant(figName);
        hgsave(h,figName);
        close(h);
    end
    
end

function [out]=OptimisationOutput_Final(paroptim,out,optimstruct)
    
    varExtract={'direction','knownOptim','objectiveName'};
    [direction,knownOptim,objectiveName]=ExtractVariables(varExtract,paroptim);
    varExtract={'axisRatio'};
    [axisRatio]=ExtractVariables(varExtract,paroptim.parametrisation);
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=out.marker;
    markerSmall=datestr(t,'_yymmddTHHMM');
    writeDirectory=[rootDir];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    CopyDiary(writeDirectory,marker)
    GenerateIterResultBinary(writeDirectory,marker,optimstruct)
    dat=GenerateOptimalSolDir(writeDirectory,markerSmall,direction,optimstruct);
    
    if strcmp(objectiveName,'CutCellFlow')
        [knownOptim]=SupersonicOptimLinRes(paroptim,rootDir,...,
            dat.xMin,dat.xMax,dat.A,dat.nPoints);
        knownOptim=0;
        tecPlotFile{1}=[writeDirectory,filesep,'Tec360plt_Flow_',marker,'.plt'];
        tecPlotFile{2}=[writeDirectory,filesep,'Tec360plt_Snak_',marker,'.plt'];
        [FID]=OpenOptimumFlowLayFile(writeDirectory,marker);
        PersnaliseLayFile(FID,tecPlotFile(2:-1:1));
        
        ExtractOptimalFlow(optimstruct,writeDirectory,direction,tecPlotFile,axisRatio,paroptim);
    end
    
    [h]=OptimHistory(optimstruct,knownOptim,direction);
    %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
    figName=[writeDirectory,'\Optimisation_',marker,'.fig'];
    figName=MakePathCompliant(figName);
    hgsave(h,figName);
   
end

function [out]=OptimisationOutput_Final_Post(paroptim,out,optimstruct)
    
    varExtract={'direction','knownOptim','objectiveName'};
    [direction,knownOptim,objectiveName]=ExtractVariables(varExtract,paroptim);
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=out.marker;
    markerSmall=datestr(t,'_yymmddTHHMM');
    writeDirectory=[rootDir];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    
    
    if strcmp(objectiveName,'CutCellFlow')

        tecPlotFile{1}=[writeDirectory,filesep,'Tec360plt_Flow_',marker,'.plt'];
        tecPlotFile{2}=[writeDirectory,filesep,'Tec360plt_Snak_',marker,'.plt'];
        [FID]=OpenOptimumFlowLayFile(writeDirectory,marker);
        PersnaliseLayFile(FID,tecPlotFile(2:-1:1));
        
        ExtractOptimalFlow(optimstruct,writeDirectory,direction,tecPlotFile,axisRatio);
    end
    
    [h]=OptimHistory(optimstruct,knownOptim,direction);
    
   
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
        str=[str,' , ', num2str(population(ii).constraint)];
        
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

%% Output tecplot video

function []=TecplotPortion_Profile(nIter,nPop,nProf,profPath,baseGrid,...
        cellCentredGrid,snaxel,snakposition)
    
    time=nIter+nProf/10^(ceil(log10(nPop)));
    fName=[profPath,filesep,'tecsubfile_',int2str(nIter),'_',int2str(nProf),'.dat'];
    FID=fopen(fName,'w');
    TecplotOutput('optim',FID,baseGrid,cellCentredGrid,snaxel,snakposition,time)
    
end

function []=TecplotPortion_Init(nIter,nPop,nProf,profPath,baseGrid,fineGrid,...
        cellCentredGrid,snaxel,snakposition)
    
    time=nIter*(nPop+1)+nProf;
    fName=[profPath,filesep,'tecsubfile_',int2str(nIter),'_',int2str(nProf),'.dat'];
    FID=fopen(fName,'w');
    TecplotOutput('optiminit',FID,baseGrid,fineGrid,cellCentredGrid,snaxel,snakposition)
    
end

function []=ConcatenateTecplotFile(iterDir,tecFilePath)
    [profPaths]=FindProfile(iterDir);
    
    compType=computer;
    
    if strcmp(compType(1:2),'PC')
        for ii=1:length(profPaths)
            system(['type "',profPaths{ii},'" >> "',tecFilePath,'"']);
        end
    else
        for ii=1:length(profPaths)
            system(['cat ''',profPaths{ii},''' >> ''',tecFilePath,'''']);
        end
    end
    
    
end

function [profPaths]=FindProfile(iterDir)
    
    [returnPath]=FindDir(iterDir,'profile',true);
    
    for ii=1:length(returnPath)
        
        profPaths(ii)=FindDir(returnPath{ii},'tecsubfile',false);
        
    end

end

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
        disp('FindDir Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
      
    
    
end

function [FID]=OpenOptimLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['tec360lay_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_ViewIterOptim.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfile(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function [FID]=OpenOptimumFlowLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['EvolOptim_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_EvolOptim.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfile(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function []=ExtractOptimalFlow(optimstruct,rootFolder,dirOptim,tecPlotFile,ratio,paramoptim)
    
    
    varExtract={'defaultVal'};
    [defaultVal]=ExtractVariables(varExtract,paramoptim);
    
    delete(tecPlotFile{1});
    delete(tecPlotFile{2});
    
    [~,iterFolders]=FindDir( rootFolder,'iteration',true);
    isIter0=regexp(iterFolders,'_0');
    
    for ii=1:length(isIter0)
        isIter0Log(ii)=~isempty(isIter0{ii});
    end
    iter0Path=[rootFolder,filesep,iterFolders{isIter0Log}];
    [initTecFolderFile]=FindDir( iter0Path,'profile',true);
    [initTecFile,initTecFileName]=FindDir(initTecFolderFile{1},'tecsubfile',false);
    jj=1;
    while ~isempty(regexp(initTecFileName{jj},'copy', 'once'))
        jj=jj+1;
    end
    [destPath]=EditVariablePLT_FEPOLYGON(1:2,[1,ratio],[0.001*pi,0.001*pi],...
        regexprep(initTecFile{jj},initTecFileName{jj},''),initTecFileName{jj},2);
    compType=computer;
    if strcmp(compType(1:2),'PC')
        system(['type "',destPath,'" > "',tecPlotFile{2},'"']);
    else
        system(['cat ''',destPath,''' > ''',tecPlotFile{2},'''']);
    end
    
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar])+defaultVal;

    for ii=1:nIter
        nAct=length(optimstruct(ii).population);
        iterRes(ii,1:nAct)=[optimstruct(ii).population(:).objective];
        
    end
    
    switch dirOptim
        case 'min'
            [minRes,minPos]=min(iterRes,[],2);
        case 'max'
            [minRes,minPos]=max(iterRes,[],2);
    end
    % Prepare CFD file with newest version
    for ii=1:nIter
        minIterPos=optimstruct(ii).population(minPos(ii)).location;
%         PrepareCFDPostProcessing(minIterPos);
        
    end
    
    kk=1;
    needRerun(kk)=1;
    for ii=2:nIter
        
        precIterPos=optimstruct(ii-1).population(minPos(ii-1)).location;
        minIterPos=optimstruct(ii).population(minPos(ii)).location;
        if ~strcmp(precIterPos,minIterPos)
            kk=kk+1;
            needRerun(kk)=ii;
            
        end
        
    end
    disp([int2str(kk), ' Reruns needed, stop bitching and be patient'])
    for jj=1:kk
        
        ii=needRerun(jj);

        minIterPos=optimstruct(ii).population(minPos(ii)).location;
        %if isempty(FindDir([minIterPos,filesep,'CFD'],'flowplt_cell',false))
        
        
            RunCFDPostProcessing(minIterPos);
            if isempty(FindDir([minIterPos,filesep,'CFD'],'flowplt_cell',false))
                CutCellFlow_Handler(paramoptim,minIterPos)
                RunCFDPostProcessing(minIterPos);
            end
        %end
        
    end
    
    for ii=1:nIter
        
        minIterPos=optimstruct(ii).population(minPos(ii)).location;
        [~,filename]=FindDir( minIterPos,'tecsubfile',false);
        jj=1;
        while ~isempty(regexp(filename{jj},'copy', 'once'))
            jj=jj+1;
        end
        
        copyfile([minIterPos,filesep,filename{jj}],[minIterPos,filesep,filename{jj},int2str(ii)])
        
        copyfile([[minIterPos,filesep,'CFD'],filesep,'flowplt_cell.plt'],...
            [[minIterPos,filesep,'CFD'],filesep,'flowplt_cell.plt',int2str(ii)])
        
        [snakPlt{ii}]=EditPLTTimeStrand(ii,3,2,minIterPos,[filename{jj},int2str(ii)]);
        [flowPlt{ii}]=EditPLTTimeStrand(ii,1,2,[minIterPos,filesep,'CFD'],...
            ['flowplt_cell.plt',int2str(ii)]);
    end
    
    for ii=1:nIter
        if strcmp(compType(1:2),'PC')
            [~,~]=system(['type "',flowPlt{ii},'" >> "',tecPlotFile{1},'"']);
            [~,~]=system(['type "',snakPlt{ii},'" >> "',tecPlotFile{2},'"']);
        else
            [~,~]=system(['cat ''',flowPlt{ii},''' >> ''',tecPlotFile{1},'''']);
            [~,~]=system(['cat ''',snakPlt{ii},''' >> ''',tecPlotFile{2},'''']);
        end
    end
    
    
    
    
    
end

function RunCFDPostProcessing(profilePath)
    
    compType=computer;
    profilePath=MakePathCompliant(profilePath);
    
    postPath=[profilePath,filesep,'CFD',filesep,'RunPost'];
    if strcmp(compType(1:2),'PC')

        [~,~]=system(['"',postPath,'.bat"']);
    else
        
        system(['echo 1 > ''',profilePath,filesep,'CFD',filesep,'nstep.txt''']);
        [~,~]=system(['''',postPath,'.sh''']);
        
    end
    
end

function PrepareCFDPostProcessing(profilePath)
    
    compType=computer;
    profilePath=MakePathCompliant(profilePath);
    
    cfdPath=[profilePath,filesep,'CFD',filesep];
    rootCode=['.',filesep,'Result_Template',filesep,'CFD_code_Template',...
        filesep,'supersonic_biplane',filesep];
    
    
    
    postPath=[profilePath,filesep,'CFD',filesep,'RunPost'];
    if strcmp(compType(1:2),'PC')
        try
            rmdir([cfdPath], 's')
        catch ME
            disp(ME.getReport)
        end
        copyfile([rootCode],[cfdPath])
        
    else
        try
            rmdir([cfdPath], 's')
        catch ME
            disp(ME.getReport)
        end
        CopyFileLinux([rootCode],[cfdPath])
        
    end
    
end

function [destPath]=EditPLTTimeStrand(time,strand,nOccur,cfdPath,filename)
    
    cfdPath=MakePathCompliant(cfdPath);
    sourceF=fopen([cfdPath,filesep,filename],'r');
    destPath=[cfdPath,filesep,'copy_',filename];
    destF=fopen(destPath,'w');
    
    flagFinished=false;
    flagTime=0;
    flagStrand=0;
    while (~feof(sourceF)) && ~flagFinished
        str=fgetl(sourceF);
        if ~isempty(regexp(str,'SOLUTIONTIME', 'once'))
            str=['SOLUTIONTIME=',num2str(time)];
            flagTime=flagTime+1;
        end
        if ~isempty(regexp(str,'STRANDID', 'once'))
            str=['STRANDID=',int2str(strand+flagStrand)];
            flagStrand=flagStrand+1;
        end
        fprintf(destF,[str,'\n']);
        flagFinished=flagTime>=nOccur && flagStrand>=nOccur;
    end
    while (~feof(sourceF)) 
        str=fgetl(sourceF);
        fprintf(destF,[str,'\n']);
    end
    fclose('all');
end

function [destPath]=EditVariablePLT_FELINESEG(varList,ratio,cfdPath,filename)
    
    cfdPath=MakePathCompliant(cfdPath);
    sourceF=fopen([cfdPath,filesep,filename],'r');
    destPath=[cfdPath,filesep,'copy_',filename];
    destF=fopen(destPath,'w');
    
    flagFinished=false;
    flagTime=0;
    flagStrand=0;
    while (~feof(sourceF)) && ~flagFinished
        str=fgetl(sourceF);
        
        fprintf(destF,[str,'\n']);
        if ~isempty(regexp(str,'NODES', 'once'))
            nNodes=regexp(str,'\d*',match);
        end
        if ~isempty(regexp(str,'FELINESEG', 'once'))
            flagFinished=true;
        end
            
    end
    for ii=1:nNodes
        str=fgetl(sourceF);
        num=str2num(str); %#ok<ST2NM>
        num(varList)=num(varList).*ratio;
        str=num2str(num);
        fprintf(destF,[str,'\n']);
    end
    while (~feof(sourceF)) 
        str=fgetl(sourceF);
        fprintf(destF,[str,'\n']);
    end
    fclose('all');
end

function [destPath]=EditVariablePLT_FEPOLYGON(varList,ratio,offsets,cfdPath,filename,nZone)
    
    cfdPath=MakePathCompliant(cfdPath);
    sourceF=fopen([cfdPath,filesep,filename],'r');
    destPath=[cfdPath,filesep,'copy_',filename];
    destF=fopen(destPath,'w');
    for jj=1:nZone
        flagFinished=false;
        flagTime=0;
        flagStrand=0;
        while (~feof(sourceF)) && ~flagFinished
            str=fgetl(sourceF);
            if ~isempty(regexp(str,'FEPOLYGON', 'once'))
                flagFinished=true;
                
            end
            fprintf(destF,[str,'\n']);
        end
        
        for ii=1:max(varList)
            str=fgetl(sourceF);
            isvar=ii==varList;
            if sum(isvar)
                
                num=str2num(str); %#ok<ST2NM>
                num=num*ratio(isvar)+offsets(isvar);
                str=num2str(num);
            end
            fprintf(destF,[str,'\n']);
        end
    end
    while (~feof(sourceF)) 
        str=fgetl(sourceF);
        fprintf(destF,[str,'\n']);
    end
    fclose('all');
end

%% 
function [h]=OptimHistory(optimstruct,knownOptim,dirOptim)
    
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    nVar=length(optimstruct(1).population);
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar]);
    hold on
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
        lSub1(1)=plot(ones(1,nVar)*ii,iterRes(ii,:),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            iterRes(iterRes==0)=1000;
            minRes=min(iterRes,[],2);
        case 'max'
            iterRes(iterRes==0)=-1000;
            minRes=max(iterRes,[],2);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:nIter,minRes,'r-');
    lSub1(3)=plot(1:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(4)=plot([0,nIter],[knownOptim knownOptim],'r--');
    
    legend(lSub1,{'Population',['Population ',dirOptim,'imum'],'Population mean',...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
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
    
    xT=box(2)-(box(2)-box(1))*0.05;
    yT=min(minRes)-(box(4)-box(3))*0.05;
    strT=['$\quad\quad$ $J^*(\mathbf{x})$ = ',sprintf('%10.3e',(min(minRes)))];
    strT={strT,['$J^*(\mathbf{x})-J^*_T$ = ',sprintf('%10.3e',min(minRes)-knownOptim)]};

    strT=regexprep(strT,'\ ','\\space');
    text(xT,yT,strT, 'interpreter','latex','HorizontalAlignment','right');
    
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
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
        lSub2(3)=semilogy([ii ii],stdline,'g:+');
    end
    legend(lSub2,{['Population ',dirOptim,'imum'],'Population mean',...
        'Standard Deviation'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J^*_T-J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    
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

function [dat]=GenerateOptimalSolDir(resultDirectory,markerSmall,optimDirection,optimstruct)
    
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
    
    dat.A=optimstruct(end).population(posOpt).additional.A;
    coord=vertcat(loop(:).subdivision);
    dat.xMin=min(coord(:,1));
    dat.xMax=max(coord(:,1));
    dat.t=max(coord(:,2))-min(coord(:,2));
    dat.nPoints=length(coord(:,1));
    
end


