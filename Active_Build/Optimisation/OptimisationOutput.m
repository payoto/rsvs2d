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

% 
function [] = OptimisationOutput()
%    FUNCTIONLIST allows local functions to be used globally once it has
%    been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory_ExecOptim'];
    HeaderActivation(funcHandles,funcDir)
    
end



function [out]=OptimisationOutput2(entryPoint,paramoptim,varargin)
    
    
    switch entryPoint
        case 'init'
            [out.marker,out.tOutput,out.rootDir]=OptimisationOutput_Init(paramoptim);
        case 'profile'
            out=varargin{1};
            out.dirprofile=OptimisationOutput_profile(paramoptim,varargin{:});
        case 'profilepost'
            out=varargin{1};
            out.dirprofile=OptimisationOutput_profilepost(paramoptim,varargin{:});
        case 'iteration'
            out=OptimisationOutput_iteration(paramoptim,varargin{:});
        case 'optstruct'
            [out]=OptimisationOutput_iterationFullPop(varargin{:},paramoptim);
        case 'final'
            out=OptimisationOutput_Final(paramoptim,varargin{:});
        case 'finalpost'
            out=OptimisationOutput_Final_Post(paramoptim,varargin{:});
    end
    
    
end


%% Entry points

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
    pltFile=['Tec360PLT_',marker,'.plt'];
    PersnaliseLayFile(FID,pltFile);
end

function [writeDirectory]=OptimisationOutput_profile(paramoptim,out,nIter,nProf,loop,...
        restartsnak,snakSave,tecStruct)
    
    varExtract={'resampleSnak','buildInternal','typeLoop'};
    [resampleSnak,buildInternal,typeLoop]=ExtractVariables(varExtract,paramoptim);
    
    
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
    if nIter==0
        savStruct.grid.base=tecStruct.baseGrid;
        savStruct.grid.refined=tecStruct.fineGrid;
        savStruct.grid.cellrefined=restartsnak.cellCentredGrid;
        savStruct.grid.connec=restartsnak.volfracconnec;
        try
            h=CheckGrid(tecStruct.baseGrid);
            hgsave(h,[writeDirectory,filesep,'gridFig.fig'])
            print(h,'-dpng','-r600',[writeDirectory,filesep,'gridFig.png'])
            close(h)
        catch ME
            disp(ME.getReport)
        end
        
    end
    
    % Output boundary data file
    [fidBoundary]=OpenBoundaryFile(writeDirectory,markerShort);
    
    if resampleSnak && ~strcmp(typeLoop,'subdivspline')
        typeLoop='subdivspline';
        warning(['Variable typeLoop was overwritten from value ''',typeLoop,''' to ''subdivspline'''])
    end
    
    
    if nIter>0
        TecplotPortion_Profile(nIter,tecStruct.nPop,nProf,writeDirectory,tecStruct.baseGrid,...
            tecStruct.volumefraction,tecStruct.snaxel,tecStruct.snakposition);
    else
        TecplotPortion_Init(nIter,tecStruct.nPop,nProf,writeDirectory,tecStruct.baseGrid,...
            tecStruct.fineGrid,tecStruct.volumefraction,tecStruct.snaxel,...
            tecStruct.snakposition);
    end
    
    GenerateProfileBinary(writeDirectory,markerShort,savStruct)
    WriteFullInfoProfile(writeDirectory,nProf,marker,t,nIter)
    
    BoundaryOutput(loop,fidBoundary,typeLoop,buildInternal);
    % DVP boundary.dat
    % or surface.xyz & displacements.xyz
    fclose(fidBoundary);
end

function [writeDirectory]=OptimisationOutput_profilepost(paramoptim,out,loop,...
        restartsnak,snakSave,tecStruct,nIter,nProf)
    
    varExtract={'resampleSnak','buildInternal','typeLoop'};
    [resampleSnak,buildInternal,typeLoop]=ExtractVariables(varExtract,paramoptim);
    
    
    marker=out.marker;
    t=out.tOutput;
    rootDir=out.rootDir;
    
    if ~(nIter==0 && nProf==1)
        iterStr=['\iteration_',int2str(nIter)];
        profStr=['\profile_',int2str(nProf)];
        writeDirectory=[rootDir,iterStr,profStr];
        markerShort=[int2str(nIter),'_',int2str(nProf)];
    else
        markerShort=['InterpOptim'];
        writeDirectory=[rootDir,filesep,'InterpOptim'];
    end
    writeDirectory=MakePathCompliant(writeDirectory);
    mkdir(writeDirectory);
    
    savStruct.restartsnak=restartsnak;
    savStruct.snakSave=snakSave;
    savStruct.loop=loop;
    
    
    % Output boundary data file
    [fidBoundary]=OpenBoundaryFile(writeDirectory,markerShort);
    
    if resampleSnak && ~strcmp(typeLoop,'subdivspline')
        typeLoop='subdivspline';
        warning(['Variable typeLoop was overwritten from value ''',typeLoop,''' to ''subdivspline'''])
    end
    
    
    
    TecplotPortion_Profile(0,1,1,writeDirectory,tecStruct.baseGrid,...
        tecStruct.volumefraction,tecStruct.snaxel,tecStruct.snakposition);
    
    
    GenerateProfileBinary(writeDirectory,markerShort,savStruct)
    WriteFullInfoProfile(writeDirectory,1,marker,t,0)
    
    BoundaryOutput(loop,fidBoundary,typeLoop,buildInternal);
    fclose(fidBoundary);
end

function [out]=OptimisationOutput_iteration(paramoptim,nIter,out,population,errorReports)
    
    t=out.tOutput;
    rootDir=out.rootDir;
    fullMark=out.marker;
    %     marker=['iteration_',int2str(nIter),datestr(t,'_yymmddTHHMM')];
    %     iterStr=['\iteration_',int2str(nIter),'_',datestr(t,'yymmddTHHMM')];
    marker=['iteration_',int2str(nIter)];
    iterStr=['\iteration_',int2str(nIter)];
    writeDirectory=[rootDir,iterStr];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    
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
        %CopyDiary(writeDirectory,marker)
        GeneratePopulationBinary(writeDirectory,marker,population)
        
        
        
        %         h=CheckOptimProfile('iter_all',writeDirectory);
        %         %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
        %         figName=[writeDirectory,'\profiles_',marker,'.fig'];
        %         figName=MakePathCompliant(figName);
        %         hgsave(h,figName);
        %         close(h);
    end
    tecFilePath=[rootDir,filesep,'Tec360PLT_',fullMark,'.plt'];
    tecFilePath=MakePathCompliant(tecFilePath);
    if ExtractVariables({'useSnake'},paramoptim)
        axisRatio=ExtractVariables({'axisRatio'},paramoptim.parametrisation);
        ConcatenateTecplotFile(writeDirectory,tecFilePath,axisRatio)
    end
end

function [out]=OptimisationOutput_iterationFullPop(out,optimstruct,paroptim)
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=out.marker;
    writeDirectory=[rootDir];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    
    GenerateIterResultBinary(writeDirectory,[marker,'partial'],optimstruct,paroptim)
    
end

function [out]=OptimisationOutput_Final(paroptim,out,optimstruct)
    
    varExtract={'direction','knownOptim','objectiveName','defaultVal',...
        'optimMethod','useSnake'};
    [direction,knownOptim,objectiveName,defaultVal,optimMethod,useSnake]...
        =ExtractVariables(varExtract,paroptim);
    varExtract={'axisRatio'};
    [axisRatio]=ExtractVariables(varExtract,paroptim.parametrisation);
    
    t=out(end).tOutput;
    rootDir=out(end).rootDir;
    allRootDir={out(:).rootDir};
    marker=out(end).marker;
    markerSmall=datestr(t,'_yymmddTHHMM');
    writeDirectory=[rootDir];
    writeDirectory=MakePathCompliant(writeDirectory);
    for ii=1:numel(optimstruct);
        for jj=1:numel(optimstruct(ii).population);
            if optimstruct(ii).population(jj).objective<=1e-16;
                optimstruct(ii).population(jj).objective=1000;
                optimstruct(ii).population(jj).constraint=0;
            end;
        end;
    end
    %CopyDiary(writeDirectory,marker)
    GenerateIterResultBinary(writeDirectory,marker,optimstruct,paroptim)
    
    
    % Figure
    [isGradient]=CheckIfGradient(optimMethod);
    isAnalytical=CheckIfAnalytical(objectiveName);
    if useSnake
        dat=GenerateOptimalSolDir(writeDirectory,markerSmall,direction,optimstruct);
        
        % Tecplot flow files
%         if strcmp(objectiveName,'CutCellFlow')
%             [knownOptim]=SupersonicOptimLinRes(paroptim,rootDir,...
%                 dat.xMin,dat.xMax,dat.A,dat.nPoints);
%         end
        
    end
    if numel(optimstruct)>=4
        [h]=OptimHistory(isGradient,optimstruct,knownOptim,defaultVal,direction);
        if ~isGradient
            [teststruct]=RebuildFromIter(optimstruct,out);
            GenerateTestResultBinary(writeDirectory,marker,teststruct,paroptim);
            [h2]=OptimHistory(isGradient,teststruct,knownOptim,defaultVal,direction);
            h=[h,h2];
        end
        
        if isAnalytical
            [h2]=PlotAnalyticalPath(paroptim,optimstruct,objectiveName);
            h=[h,h2];
        end
        %print(h,'-depsc','-r600',[writeDirectory,'\profiles_',marker,'.eps']);
        if ~isGradient
            figList={'\Optimisation_','\DesVarHist_','\OptimisationTest_','\DesVarHistTest_'};
        else
            figList={'\Optimisation_','\GradHistory_'};
            [h2]=PlotGradientsRefine(paroptim,optimstruct);
            h=[h,h2];
        end
        if isAnalytical
            figList{end+1}={'\DesVarPath'};
        else
            %figList{end+1}={'\DesVarPath'};
        end
        
        for ii=1:length(h)
            figName=[writeDirectory,figList{ii},marker,'.fig']
            figName=MakePathCompliant(figName);
            hgsave(h(ii),figName);
        end
    end
    if strcmp(objectiveName,'CutCellFlow')

        tecPlotFile{1}=['Tec360plt_Flow_',marker,'.plt'];
        tecPlotFile{2}=['Tec360plt_Snak_',marker,'.plt'];
        tecPlotPre{1}=['Tec360plt_Flow_',marker,'_pre.plt'];
        tecPlotPre{2}=['Tec360plt_Snak_',marker,'_pre.plt'];
        [FID]=OpenOptimumFlowLayFile(writeDirectory,marker);
        PersnaliseLayFile(FID,tecPlotPre(2:-1:1));
        tecPlotFile{1}=[writeDirectory,filesep,tecPlotFile{1}];
        tecPlotFile{2}=[writeDirectory,filesep,tecPlotFile{2}];
        ExtractOptimalFlow(optimstruct,writeDirectory,direction,...
            tecPlotFile,axisRatio,paroptim,allRootDir,isGradient);
        
    end
    
    if strcmp(objectiveName,'ASOFlow')

        tecPlotFile{1}=['Tec360plt_Flow_',marker,'.plt'];
        tecPlotFile{2}=['Tec360plt_Snak_',marker,'.plt'];
        tecPlotPre{1}=['Tec360plt_Flow_',marker,'_pre.plt'];
        tecPlotPre{2}=['Tec360plt_Snak_',marker,'_pre.plt'];
        [FID]=OpenOptimumFlowLayFile(writeDirectory,marker);
        PersnaliseLayFile(FID,tecPlotPre(2:-1:1));
        tecPlotFile{1}=[writeDirectory,filesep,tecPlotFile{1}];
        tecPlotFile{2}=[writeDirectory,filesep,tecPlotFile{2}];
        ExtractOptimalFlowSU2(optimstruct,writeDirectory,direction,...
            tecPlotFile,axisRatio,paroptim,allRootDir,isGradient);
        tecPlotFile{1}=regexprep(tecPlotFile{1},'Flow','FlowMulti');
        tecPlotFile{2}=regexprep(tecPlotFile{2},'Snak','SnakMulti');
        ExtractOptimalFlowSU2Multi(optimstruct,writeDirectory,direction,...
            tecPlotFile,axisRatio,paroptim,allRootDir,isGradient,...
            max(cellfun(@numel,{optimstruct.population})));
        
    end
    
    if strcmp(objectiveName,'InverseDesign') || strcmp(objectiveName,'InverseDesignTopo')
        tecPlotFile{1}=['Tec360plt_Flow_',marker,'.plt'];
        tecPlotFile{2}=['Tec360plt_Snak_',marker,'.plt'];
        tecPlotPre{1}=['Tec360plt_Flow_',marker,'_pre.plt'];
        tecPlotPre{2}=['Tec360plt_Snak_',marker,'_pre.plt'];
        
        [FID]=OpenOptimumSnakLayFile(writeDirectory,marker);
        PersnaliseLayFile(FID,tecPlotPre(2:-1:1));
        tecPlotFile{1}=[writeDirectory,filesep,tecPlotFile{1}];
        tecPlotFile{2}=[writeDirectory,filesep,tecPlotFile{2}];
        ExtractOptimalSnake(optimstruct,writeDirectory,direction,...
            tecPlotFile,axisRatio,paroptim,allRootDir);
        
    end
end

function [out]=OptimisationOutput_Final_Post(paroptim,out,optimstruct)
    
    
    error('Function Is Deprecated')
    
    varExtract={'direction','knownOptim','objectiveName'};
    [direction,knownOptim,objectiveName]=ExtractVariables(varExtract,paroptim);
    
    t=out.tOutput;
    rootDir=out.rootDir;
    marker=out.marker;
    markerSmall=datestr(t,'_yymmddTHHMM');
    writeDirectory=[rootDir];
    writeDirectory=MakePathCompliant(writeDirectory);
    
    
    
    if strcmp(objectiveName,'CutCellFlow')
        
        tecPlotFile{1}=['Tec360plt_Flow_',marker,'.plt'];
        tecPlotFile{2}=['Tec360plt_Snak_',marker,'.plt'];
        [FID]=OpenOptimumFlowLayFile(writeDirectory,marker);
        PersnaliseLayFile(FID,tecPlotFile(2:-1:1));
        tecPlotFile{1}=[writeDirectory,filesep,tecPlotFile{1}];
        tecPlotFile{1}=[writeDirectory,filesep,tecPlotFile{2}];
        
        ExtractOptimalFlow(optimstruct,writeDirectory,direction,tecPlotFile,axisRatio);
    end
    
    [h]=OptimHistory(optimstruct,knownOptim(1),direction);
    
    
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
    entryCell{1}=[entryCell{1},', exception , fill , nonfillvar'];
    
    for ii=1:length(population)
        
        str=[int2str(nIter)];
        str=[str,' , ', int2str(ii)];
        str=[str,' , ', num2str(population(ii).objective,' %30.24e ')];
        str=[str,' , ', num2str(population(ii).constraint)];
        
        for jj=1:nAdditional
            try 
                str=[str,' , ', num2str(population(ii).additional.(fieldsAdd{jj}),' %30.24e ')];
            catch
                str=[str,' , ',''];
            end
        end
        str=[str,' , ', population(ii).exception];
        str=[str,' , ', num2str(population(ii).fill,' %30.24e ')];
        str=[str,' , ', num2str(population(ii).nonfillvar,' %30.24e ')];
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

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    %     if iscell(rootDir)
    %         subDir=dir(rootDir{1});
    %         subDir(1:2)=[];
    %         for ii=2:numel(rootDir)
    %             partsubDir=dir(rootDir{ii});
    %             partsubDir(1:2)=[];
    %             subDir=[subDir,partsubDir];
    %         end
    %     else
    subDir=dir(rootDir);
    if isempty(subDir)
        errStruct.message=sprintf('rootDir is an empty directory structure\n %s',rootDir);
        errStruct.identifier='Optimiser:Output:findDirFail';
        error(errStruct)
    else
        subDir(1:2)=[];
        %     end
        nameCell={subDir(:).name};
        isprofileDirCell=strfind(nameCell,strDir);
        for ii=1:length(subDir)
            subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
                ~xor(subDir(ii).isdir,isTargDir);
        end
        
        returnSub=find([subDir(:).isProfile]);
        
        
        if isempty(returnSub)
            fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
        end
        for ii=1:length(returnSub)
            returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
            returnName{ii}=subDir(returnSub(ii)).name;
            
        end
    end
    
    
end

function [FID]=OpenOptimLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['tec360lay_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_ViewIterOptim.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfileRobust(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function [FID]=OpenOptimumFlowLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['EvolOptim_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_EvolOptim.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfileRobust(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function [FID]=OpenOptimumSnakLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['EvolSnake_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_OptimSnak.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfileRobust(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function [rootDirName]=InitOptimalFlowOutput(rootFolder,ratio,tecPlotFile)
    try
        [~,iterFolders]=FindDir( rootFolder,'iteration',true);
    catch
        warning('%s does not exist anymore',rootFolder)
        rootDirName='';
        iterFolders={};
    end
    
    if ~isempty(iterFolders)
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
        [destPath]=EditVariablePLT_FEPOLYGON(1:2,[1,ratio],[1e-6*pi/3.14,1e-6*pi/3.14],...
            regexprep(initTecFile{jj},initTecFileName{jj},''),initTecFileName{jj},2);

        compType=computer;
        if strcmp(compType(1:2),'PC')
            system(['type "',destPath,'" >> "',tecPlotFile{2},'"']);
        else
            system(['cat ''',destPath,''' >> ''',tecPlotFile{2},'''']);
        end
        listFileSep=regexp(rootFolder,filesep);
        rootDirName=rootFolder(listFileSep(end)+1:end);
    end
end

function [tecPlotPre]=ExtractOptimalFlow(optimstruct,rootFolder,dirOptim,...
        tecPlotFile,ratio,paramoptim,allRootDir,isGradient)
    
    
    varExtract={'defaultVal','worker','CFDfolder'};
    [defaultVal,worker,CFDfolder]=ExtractVariables(varExtract,paramoptim);
    
    delete(tecPlotFile{1});
    delete(tecPlotFile{2});
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    
    for ii=1:numel(allRootDir)
        rootDirName{ii}=InitOptimalFlowOutput(allRootDir{ii},ratio,tecPlotFile);
    end
    rootDirName=rootDirName(~cellfun(@isempty,rootDirName));
    compType=computer;
    for ii=1:nIter

        itL=[optimstruct(ii).population(:).objective];
        nVarLoc=length(itL);
        iterRes(ii,1:nVarLoc)=itL;
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);

    end
    
    switch dirOptim
        case 'min'
            [minRes,minPos]=min(iterRes,[],2);
        case 'max'
            [minRes,minPos]=max(iterRes,[],2);

    end
    if isGradient
        minPos(1:2:end)=1;
        minRes(1:2:end)=iterRes(1:2:end,1);
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
    %parfor jj=1:kk
    postList=needRerun;
    postLog=true(size(1:kk));
    parfor jj=1:kk
        ii=needRerun(jj);
        minIterPos=optimstruct(ii).population(minPos(ii)).location;

        try

            if isempty(FindDir([minIterPos,filesep,'CFD'],'flowplt_cell',false))
                RunCFDPostProcessing(minIterPos);
                if isempty(FindDir([minIterPos,filesep,'CFD'],'flowplt_cell',false))
                    CutCellFlow_Handler(paramoptim,minIterPos)
                    RunCFDPostProcessing(minIterPos);
                    if isempty(FindDir([minIterPos,filesep,'CFD'],'flowplt_cell',false))
                        PrepareCFDPostProcessing(minIterPos,CFDfolder);
                        CutCellFlow_Handler(paramoptim,minIterPos)
                        RunCFDPostProcessing(minIterPos);
                        if isempty(FindDir([minIterPos,filesep,'CFD'],'flowplt_cell',false))
                            error(sprintf('Could Not gernerate flow plot for %s',minIterPos))
                        end
                    end
                end
            end

        catch ME
            disp(ME.getReport);
            postLog(jj)=false;
        end
        
    end
    postList((~postLog))=[];
    minIterRootDirNum=zeros([1,nIter]);
    
    for ii=postList
        
        minIterPos=optimstruct(ii).population(minPos(ii)).location;
        [~,filename]=FindDir( minIterPos,'tecsubfile',false);
        jj=1;
        while isempty(regexp(minIterPos,rootDirName{jj}, 'once')) && jj<numel(rootDirName)
            jj=jj+1;
        end
        minIterRootDirNum(ii)=jj;
        jj=1;
        try
            while ~isempty(regexp(filename{jj},'copy', 'once'))
                jj=jj+1;
            end
        catch ME
            ME.getReport
            minIterPos
            optimstruct(ii).population
            jj=jj
            ii=ii
            filename
            throw(ME)
        end
        
        
        copyfileRobust([minIterPos,filesep,filename{jj}],[minIterPos,filesep,filename{jj},int2str(ii)])
        
        copyfileRobust([[minIterPos,filesep,'CFD'],filesep,'flowplt_cell.plt'],...
            [[minIterPos,filesep,'CFD'],filesep,'flowplt_cell.plt',int2str(ii)])
        
        %[snakPlt{ii}]=EditPLTTimeStrand(ii,3,2,minIterPos,[filename{jj},int2str(ii)]);
        dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'};
        expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'};
        val={ii,[3 4],minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5};
        nOccur=[2 2 1 1];
        [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur);
        
        [flowPlt{ii}]=EditPLTTimeStrand(ii,1,2,[minIterPos,filesep,'CFD'],...
            ['flowplt_cell.plt',int2str(ii)]);
    end
    
    % dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'}
    % expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'}
    % val={ii,3,minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5}
    % nOccur=[2 2 1 1]
    % [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur)
    
    for ii=postList
        
        if strcmp(compType(1:2),'PC')
            [~,~]=system(['type "',flowPlt{ii},'" >> "',tecPlotFile{1},'"']);
            [~,~]=system(['type "',snakPlt{ii},'" >> "',tecPlotFile{2},'"']);
        else
            [~,~]=system(['cat ''',flowPlt{ii},''' >> ''',tecPlotFile{1},'''']);
            [~,~]=system(['cat ''',snakPlt{ii},''' >> ''',tecPlotFile{2},'''']);
        end
    end
    tecPlotPre=regexprep(tecPlotFile,'\.plt','_pre.plt');
    if strcmp(compType(1:2),'PC')
        [d,c]=system(['preplot "',tecPlotFile{1},'" "',tecPlotPre{1},'"']);
        [d,c]=system(['preplot "',tecPlotFile{2},'" "',tecPlotPre{2},'"']);
    end
    
    
    
    
end

function [tecPlotPre]=ExtractOptimalFlowSU2(optimstruct,rootFolder,dirOptim,...
        tecPlotFile,ratio,paramoptim,allRootDir,isGradient)
    
    
    varExtract={'defaultVal','worker','CFDfolder'};
    [defaultVal,worker,CFDfolder]=ExtractVariables(varExtract,paramoptim);
    
    delete(tecPlotFile{1});
    delete(tecPlotFile{2});
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    
    for ii=1:numel(allRootDir)
        rootDirName{ii}=InitOptimalFlowOutput(allRootDir{ii},ratio,tecPlotFile);
    end
    rootDirName=rootDirName(~cellfun(@isempty,rootDirName));
    compType=computer;
    for ii=1:nIter

        itL=[optimstruct(ii).population(:).objective];
        nVarLoc=length(itL);
        iterRes(ii,1:nVarLoc)=itL;
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);

    end
    
    switch dirOptim
        case 'min'
            [minRes,minPos]=min(iterRes,[],2);
        case 'max'
            [minRes,minPos]=max(iterRes,[],2);

    end
    if isGradient
        minPos(1:2:end)=1;
        minRes(1:2:end)=iterRes(1:2:end,1);
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
    %parfor jj=1:kk
    postList=needRerun;
    postLog=true(size(1:kk));
    for jj=1:kk
        ii=needRerun(jj);
        minIterPos=optimstruct(ii).population(minPos(ii)).location;

        try
            
            fileOrig=[minIterPos,filesep,'run',filesep,'flow.dat'];
            if ~exist(fileOrig,'file')
                configOrig=[minIterPos,filesep,'run',filesep,'su2.cfg'];
                configEdit=[minIterPos,filesep,'run',filesep,'su2pp.cfg'];
                fidOrig=fopen(configOrig,'r');
                fidEdit=fopen(configEdit,'w');
                while ~feof(fidOrig)
                    fprintf(fidEdit,'%s\n',regexprep(fgetl(fidOrig),...
                        '^MATH_PROBLEM=.*$','MATH_PROBLEM=DIRECT'));
                end
                fclose(fidOrig);
                fclose(fidEdit);
                system(['cd "',minIterPos,'" && SU2_SOL run',filesep,'su2pp.cfg']);
                if ~exist(fileOrig,'file')
                    error('Could not generate solution file')
                end
            end
            if isempty(FindDir([minIterPos,filesep,'run'],'flowplt_cell',false))
                
                fileNew=[minIterPos,filesep,'run',filesep,'flowplt_cell.plt'];
                fidOrig=fopen(fileOrig,'r');
                fidNew=fopen(fileNew,'w');
                str=fgetl(fidOrig);
                fprintf(fidNew,'%s \n',fgetl(fidOrig));
                fprintf(fidNew,'ZONE \n');
                fprintf(fidNew,'STRANDID = 1 \n');
                fprintf(fidNew,'SOLUTIONTIME = 1 \n');
                fprintf(fidNew,'%s \n',regexprep(fgetl(fidOrig),'ZONE ',''));
                while ~feof(fidOrig)
                    fprintf(fidNew,'%s \n',fgetl(fidOrig));
                end
                fclose(fidOrig);
                fclose(fidNew);
            end
            
        catch ME
            disp(ME.getReport);
            postLog(jj)=false;
        end
        
    end
    postList((~postLog))=[];
    minIterRootDirNum=zeros([1,nIter]);
    
    for ii=postList
        
        minIterPos=optimstruct(ii).population(minPos(ii)).location;
        [~,filename]=FindDir( minIterPos,'tecsubfile',false);
        jj=1;
        while isempty(regexp(minIterPos,rootDirName{jj}, 'once')) && jj<numel(rootDirName)
            jj=jj+1;
        end
        minIterRootDirNum(ii)=jj;
        jj=1;
        try
            while ~isempty(regexp(filename{jj},'copy', 'once'))
                jj=jj+1;
            end
        catch ME
            ME.getReport
            minIterPos
            optimstruct(ii).population
            jj=jj
            ii=ii
            filename
            throw(ME)
        end
        
        
        copyfileRobust([minIterPos,filesep,filename{jj}],[minIterPos,filesep,filename{jj},int2str(ii)])
        
        copyfileRobust([[minIterPos,filesep,'run'],filesep,'flowplt_cell.plt'],...
            [[minIterPos,filesep,'run'],filesep,'flowplt_cell.plt',int2str(ii)])
        
        %[snakPlt{ii}]=EditPLTTimeStrand(ii,3,2,minIterPos,[filename{jj},int2str(ii)]);
        dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'};
        expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'};
        val={ii,[3 4],minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5};
        nOccur=[2 2 1 1];
        [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur);
        
        [flowPlt{ii}]=EditPLTTimeStrand(ii,1,2,[minIterPos,filesep,'run'],...
            ['flowplt_cell.plt',int2str(ii)]);
    end
    
    % dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'}
    % expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'}
    % val={ii,3,minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5}
    % nOccur=[2 2 1 1]
    % [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur)
    
    for ii=postList
        
        if strcmp(compType(1:2),'PC')
            [~,~]=system(['type "',flowPlt{ii},'" >> "',tecPlotFile{1},'"']);
            [~,~]=system(['type "',snakPlt{ii},'" >> "',tecPlotFile{2},'"']);
        else
            [~,~]=system(['cat ''',flowPlt{ii},''' >> ''',tecPlotFile{1},'''']);
            [~,~]=system(['cat ''',snakPlt{ii},''' >> ''',tecPlotFile{2},'''']);
        end
    end
    tecPlotPre=regexprep(tecPlotFile,'\.plt','_pre.plt');
    if strcmp(compType(1:2),'PC')
        [d,c]=system(['preplot "',tecPlotFile{1},'" "',tecPlotPre{1},'"']);
        [d,c]=system(['preplot "',tecPlotFile{2},'" "',tecPlotPre{2},'"']);
    end
    
    
    
    
end


function [tecPlotPre]=ExtractOptimalFlowSU2Multi(optimstruct,rootFolder,dirOptim,...
        tecPlotFile,ratio,paramoptim,allRootDir,isGradient,n)
    
    
    varExtract={'defaultVal','worker','CFDfolder'};
    [defaultVal,worker,CFDfolder]=ExtractVariables(varExtract,paramoptim);
    
    delete(tecPlotFile{1});
    delete(tecPlotFile{2});
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    
    for ii=1:numel(allRootDir)
        rootDirName{ii}=InitOptimalFlowOutput(allRootDir{ii},ratio,tecPlotFile);
    end
    rootDirName=rootDirName(~cellfun(@isempty,rootDirName));
    compType=computer;
    switch dirOptim
        case 'min'
            dir='ascend';
        case 'max'
            dir='descend';

    end
    minPos=zeros([2 0]);
    for ii=1:nIter

        itL=[optimstruct(ii).population(:).objective];
        [~,itOrd]=sort(itL,dir);
        nVarLoc=length(itL);
        iterRes(ii,1:nVarLoc)=itL;
        minPos=[minPos,[ii*ones([1 min(n,nVarLoc)]);...
            reshape(itOrd(1:min(n,nVarLoc)),[1 min(n,nVarLoc)])]];
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);

    end
    
    % Prepare CFD file with newest version
    for ii=1:size(minPos,2)
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;
        %         PrepareCFDPostProcessing(minIterPos);
        
    end
    
    kk=1;
    needRerun(kk)=1;
    for ii=2:size(minPos,2)
        
        precIterPos=optimstruct(minPos(1,ii-1)).population(minPos(2,ii-1)).location;
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;
        if ~strcmp(precIterPos,minIterPos)
            kk=kk+1;
            needRerun(kk)=ii;
            
        end
        
    end
    disp([int2str(kk), ' Reruns needed, stop bitching and be patient'])
    %parfor jj=1:kk
    postList=needRerun;
    postLog=true(size(1:kk));
    for jj=1:kk
        ii=needRerun(jj);
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;

        try
            
            fileOrig=[minIterPos,filesep,'run',filesep,'flow.dat'];
            if ~exist(fileOrig,'file')
                configOrig=[minIterPos,filesep,'run',filesep,'su2.cfg'];
                configEdit=[minIterPos,filesep,'run',filesep,'su2pp.cfg'];
                fidOrig=fopen(configOrig,'r');
                fidEdit=fopen(configEdit,'w');
                while ~feof(fidOrig)
                    fprintf(fidEdit,'%s\n',regexprep(fgetl(fidOrig),...
                        '^MATH_PROBLEM=.*$','MATH_PROBLEM=DIRECT'));
                end
                fclose(fidOrig);
                fclose(fidEdit);
                system(['cd "',minIterPos,'" && SU2_SOL run',filesep,'su2pp.cfg']);
                if ~exist(fileOrig,'file')
                    error('Could not generate solution file')
                end
            end
            if isempty(FindDir([minIterPos,filesep,'run'],'flowplt_cell',false))
                
                fileNew=[minIterPos,filesep,'run',filesep,'flowplt_cell.plt'];
                fidOrig=fopen(fileOrig,'r');
                fidNew=fopen(fileNew,'w');
                str=fgetl(fidOrig);
                fprintf(fidNew,'%s \n',fgetl(fidOrig));
                fprintf(fidNew,'ZONE \n');
                fprintf(fidNew,'STRANDID = 1 \n');
                fprintf(fidNew,'SOLUTIONTIME = 1 \n');
                fprintf(fidNew,'%s \n',regexprep(fgetl(fidOrig),'ZONE ',''));
                while ~feof(fidOrig)
                    fprintf(fidNew,'%s \n',fgetl(fidOrig));
                end
                fclose(fidOrig);
                fclose(fidNew);
            end
            
        catch ME
            disp(ME.getReport);
            postLog(jj)=false;
        end
        
    end
    postList((~postLog))=[];
    minIterRootDirNum=zeros([1,size(minPos,2)]);
    
    for ii=postList
        
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;
        [~,filename]=FindDir( minIterPos,'tecsubfile',false);
        jj=1;
        while isempty(regexp(minIterPos,rootDirName{jj}, 'once')) && jj<numel(rootDirName)
            jj=jj+1;
        end
        minIterRootDirNum(ii)=jj;
        jj=1;
        try
            while ~isempty(regexp(filename{jj},'copy', 'once'))
                jj=jj+1;
            end
        catch ME
            ME.getReport
            minIterPos
            optimstruct(minPos(1,ii)).population
            jj=jj
            ii=ii
            filename
            throw(ME)
        end
        
        
        copyfileRobust([minIterPos,filesep,filename{jj}],[minIterPos,filesep,...
            filename{jj},num2str(minPos(:,ii)','%i_')])
        
        copyfileRobust([[minIterPos,filesep,'run'],filesep,'flowplt_cell.plt'],...
            [[minIterPos,filesep,'run'],filesep,'flowplt_cell.plt',num2str(minPos(:,ii)','%i_')])
        
        %[snakPlt{ii}]=EditPLTTimeStrand(ii,3,2,minIterPos,[filename{jj},int2str(ii)]);
        
        decNumIterProf=(minPos(1,ii)+((minPos(2,ii))/(10^(...
            1+ceil(log10(n))+1+floor(log10(numel(optimstruct(minPos(1,ii)).population))))))...
            +(mod(ii-1,n)+1)/10^(1+ceil(log10(n))))*(10^(0*...
            max(-4+1+ceil(log10(n))+1+floor(log10(numel(optimstruct(minPos(1,ii)).population))),0)));
        dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'};
        expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'};
        val={decNumIterProf,[3 4],minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5};
        nOccur=[2 2 1 1];
        [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},num2str(minPos(:,ii)','%i_')],dat,expr,val,nOccur);
        
        [flowPlt{ii}]=EditPLTTimeStrand(decNumIterProf,1,2,[minIterPos,filesep,'run'],...
            ['flowplt_cell.plt',num2str(minPos(:,ii)','%i_')]);
    end
    
    % dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'}
    % expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'}
    % val={ii,3,minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5}
    % nOccur=[2 2 1 1]
    % [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur)
    
    for ii=postList
        
        if strcmp(compType(1:2),'PC')
            [~,~]=system(['type "',flowPlt{ii},'" >> "',tecPlotFile{1},'"']);
            [~,~]=system(['type "',snakPlt{ii},'" >> "',tecPlotFile{2},'"']);
        else
            [~,~]=system(['cat ''',flowPlt{ii},''' >> ''',tecPlotFile{1},'''']);
            [~,~]=system(['cat ''',snakPlt{ii},''' >> ''',tecPlotFile{2},'''']);
        end
    end
    tecPlotPre=regexprep(tecPlotFile,'\.plt','_pre.plt');
    if strcmp(compType(1:2),'PC')
        [d,c]=system(['preplot "',tecPlotFile{1},'" "',tecPlotPre{1},'"']);
        [d,c]=system(['preplot "',tecPlotFile{2},'" "',tecPlotPre{2},'"']);
    end
    
    
    
    
end

function [tecPlotPre]=ExtractOptimalSnake(optimstruct,rootFolder,dirOptim,...
        tecPlotFile,ratio,paramoptim,allRootDir)
    
    
    varExtract={'defaultVal','worker'};
    [defaultVal,worker]=ExtractVariables(varExtract,paramoptim);
    
    delete(tecPlotFile{1});
    delete(tecPlotFile{2});
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    
    for ii=1:numel(allRootDir)
        rootDirName{ii}=InitOptimalFlowOutput(allRootDir{ii},ratio,tecPlotFile);
    end
    rootDirName=rootDirName(~cellfun(@isempty,rootDirName));
    
    compType=computer;
    for ii=1:nIter
        itL=[optimstruct(ii).population(:).objective];
        nVarLoc=length(itL);
        iterRes(ii,1:nVarLoc)=itL;
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
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
    
    minIterRootDirNum=zeros([1,nIter]);
    
    for ii=1:nIter
        
        minIterPos=optimstruct(ii).population(minPos(ii)).location;
        [~,filename]=FindDir( minIterPos,'tecsubfile',false);
        jj=1;
        while isempty(regexp(minIterPos,rootDirName{jj}, 'once')) && jj<numel(rootDirName)
            jj=jj+1;
        end
        minIterRootDirNum(ii)=jj;
        jj=1;
        try
            while ~isempty(regexp(filename{jj},'copy', 'once'))
                jj=jj+1;
            end
        catch ME
            ME.getReport
            minIterPos
            optimstruct(ii).population
            jj=jj
            ii=ii
            filename
            throw(ME)
        end
        
        
        copyfileRobust([minIterPos,filesep,filename{jj}],[minIterPos,filesep,filename{jj},int2str(ii)])
        
%         copyfileRobust([[minIterPos,filesep,'CFD'],filesep,'flowplt_cell.plt'],...
%             [[minIterPos,filesep,'CFD'],filesep,'flowplt_cell.plt',int2str(ii)])
        
        %[snakPlt{ii}]=EditPLTTimeStrand(ii,3,2,minIterPos,[filename{jj},int2str(ii)]);
        dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'};
        expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'};
        val={ii,[3 4],minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5};
        nOccur=[2 2 1 1];
        [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur);
        
%         [flowPlt{ii}]=EditPLTTimeStrand(ii,1,2,[minIterPos,filesep,'CFD'],...
%             ['flowplt_cell.plt',int2str(ii)]);
    end
    
    % dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'}
    % expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'}
    % val={ii,3,minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5}
    % nOccur=[2 2 1 1]
    % [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur)
    
    for ii=1:nIter
        if strcmp(compType(1:2),'PC')
            %[~,~]=system(['type "',flowPlt{ii},'" >> "',tecPlotFile{1},'"']);
            [~,~]=system(['type "',snakPlt{ii},'" >> "',tecPlotFile{2},'"']);
        else
            %[~,~]=system(['cat ''',flowPlt{ii},''' >> ''',tecPlotFile{1},'''']);
            [~,~]=system(['cat ''',snakPlt{ii},''' >> ''',tecPlotFile{2},'''']);
        end
    end
    tecPlotPre=regexprep(tecPlotFile,'\.plt','_pre.plt');
    if strcmp(compType(1:2),'PC')
        %[d,c]=system(['preplot "',tecPlotFile{1},'" "',tecPlotPre{1},'"']);
        [d,c]=system(['preplot "',tecPlotFile{2},'" "',tecPlotPre{2},'"']);
    end
    
    
    
    
end

function [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal)
    
    nVar=0;
    for ii=1:length(optimstruct)
        nVar=max([length(optimstruct(ii).population),nVar]);
    end
    nIter=length(optimstruct);
    iterRes=zeros([nIter,nVar])+sign(defaultVal)*abs(defaultVal^2);
    
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

function PrepareCFDPostProcessing(profilePath,CFDfolder)
    
    compType=computer;
    profilePath=MakePathCompliant(profilePath);
    
    cfdPath=[profilePath,filesep,'CFD',filesep];
    rootCode=MakePathCompliant(CFDfolder);
    
    
    
    postPath=[profilePath,filesep,'CFD',filesep,'RunPost'];
    if strcmp(compType(1:2),'PC')
        try
            rmdir([cfdPath], 's')
        catch ME
            disp(ME.getReport)
        end
        copyfileRobust([rootCode],[cfdPath])
        
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
            str=['SOLUTIONTIME=',num2str(time,'%.8f')];
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

function [destPath]=EditPLTHeader(cfdPath,filename,dat,expr,val,nOccur)
    
    % dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'}
    % expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'}
    % val={time,strand,conneczone,conneczone}
    % nOccur=[2 2 1 1]
    %
    
    cfdPath=MakePathCompliant(cfdPath);
    sourceF=fopen([cfdPath,filesep,filename],'r');
    destPath=[cfdPath,filesep,'copy_',filename];
    destF=fopen(destPath,'w');
    
    flagFinished=false;
    jj=ones(size(dat));
    nOccurCurr=zeros(size(dat));
    while (~feof(sourceF)) && ~flagFinished
        str=fgetl(sourceF);
        
        for ii=1:numel(dat)
            if ~isempty(regexp(str,dat{ii},'once'))
                nOccurCurr(ii)=nOccurCurr(ii)+1;
                str=sprintf(expr{ii},val{ii}(min(nOccurCurr(ii),numel(val{ii}))));
                break
            end
        end
        
        fprintf(destF,[str,'\n']);
        flagFinished= all(nOccurCurr>=nOccur);
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

%% Figures

function [h]=OptimHistory(isGradient,optimstruct,knownOptim,defaultVal,dirOptim)
    
    if ~isGradient
        
        [h]=OptimHistory_nograd(optimstruct,knownOptim,defaultVal,dirOptim);
    else
        [h]=OptimHistory_grad(optimstruct,knownOptim,defaultVal,dirOptim);
    end
    
end

function [h]=OptimHistory_old(optimstruct,knownOptim,dirOptim)
    
    
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
        lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
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
    yT=(min(minRes)+max(minRes))/2;%min(minRes)-(box(4)-box(3))*0.05;
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

function [h]=OptimHistory_nograd(optimstruct,knownOptim,defaultVal,dirOptim)
    
    knownOptim=knownOptim(1);
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    hold on
    
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
        lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            [minRes,minPos]=min(iterRes,[],2);
        case 'max'
            [minRes,minPos]=max(iterRes,[],2);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:nIter,minRes,'r-');
    lSub1(3)=plot(1:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(4)=plot([0,nIter],[knownOptim(1) knownOptim(1)],'r--');
    
    legend(lSub1,{'Population',['Population ',dirOptim,'imum'],'Population mean',...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(meanRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim(1)]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(meanRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim(1)]);
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
    
    % Plot 2
    axh=subplot(1,2,2,'ticklabelinterpreter','latex');
    
    switch dirOptim
        case 'min'
            errMeasure=-(knownOptim-minRes);
            
        case 'max'
            errMeasure=knownOptim-minRes;
    end
    
    
    popNominal=zeros([nIter,numel(optimstruct(1).population(1).fill)]);
    objNominal=zeros([nIter,1]);
    cBounds=[min(minRes),max(minRes)];
    [datCol]=ProjectColormap(h.Colormap,meanRes,cBounds);
    maxFill=popNominal;
    minFill=popNominal;
    meanFill=popNominal;
    stdFill=popNominal;
    %             plot3(ones(size(optimstruct(ii).population(jj).fill))*((ii+1)/2),...
    %                 1:numel(optimstruct(ii).population(jj).fill),...
    %                 optimstruct(ii).population(:).fill,'Color',datCol(ii,:));
    for ii=1:nIter
        
        popNominal(ii,:)=optimstruct(ii).population(minPos(ii)).fill;
        objNominal(ii)=optimstruct(ii).population(minPos(ii)).objective;
        fillPop=vertcat(optimstruct(ii).population(:).fill);
        maxFill(ii,:)=max(fillPop);
        minFill(ii,:)=min(fillPop);
        meanFill(ii,:)=mean(fillPop);
        stdFill(ii,:)=std(fillPop);
        
        plot3(ones(size(optimstruct(ii).population(1).fill))*((ii+1)),...
            1:numel(optimstruct(ii).population(1).fill),...
            maxFill(ii,:),'Color',datCol(ii,:));
        hold on
        plot3(ones(size(optimstruct(ii).population(1).fill))*((ii+1)),...
            1:numel(optimstruct(ii).population(1).fill),...
            minFill(ii,:),'Color',datCol(ii,:));
    end
    
    
    iterVec=1:nIter;
    varVec=1:length(popNominal(1,:));
    
    [iterGrid,varGrid]=meshgrid(iterVec,varVec);
    [objGrid,~]=meshgrid(objNominal,varVec);
    [meanGrid,~]=meshgrid(meanRes,varVec);
    
    surf(iterGrid,varGrid,popNominal',objGrid)
    colorbar;
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    zlabel('Variable Value', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    %     axh=subplot(2,2,4,'ticklabelinterpreter','latex');
    %     hold on
    %
    %     mesh(iterGrid,varGrid,maxFill')
    %     hold on
    %     surf(iterGrid,varGrid,minFill')
    
    % Figure 2
    h(2)=figure('Name','Design Variable Evolution','Position',[20 100 1000 600]);
    
    axh=subplot(2,2,1,'ticklabelinterpreter','latex');
    s(1)=surf(iterGrid,varGrid,meanFill');
    caxis([0 1])
    
    axc(1)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    set(axh,'ticklabelinterpreter','latex')
    axh=subplot(2,2,3,'ticklabelinterpreter','latex');
    s(2)=surf(iterGrid,varGrid,stdFill');
    caxis([0 1])
    axc(2)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    set(axh,'ticklabelinterpreter','latex')
    axh=subplot(2,2,2,'ticklabelinterpreter','latex');
    s(3)=surf(iterGrid,varGrid,maxFill');
    caxis([0 1])
    axc(3)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    axh=subplot(2,2,4,'ticklabelinterpreter','latex');
    s(4)=surf(iterGrid,varGrid,minFill');
    caxis([0 1])
    axc(4)=colorbar;
    view(0 ,90)
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    
    cLabel={'Mean','Standard Deviation','Maximum','Minimum'};
    for ii=1:length(s)
        s(ii).EdgeColor='none';
        axc(ii).TickLabelInterpreter='latex';
        axc(ii).Label.String=cLabel{ii};
        axc(ii).Label.FontSize=14;
    end
    
    
    
    %{
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
    %}
    
    
end

function [h]=OptimHistory_grad(optimstruct,knownOptim,defaultVal,dirOptim)
    
    h=figure('Name','Optimisation Result','Position',[20 100 1000 600]);
    
    % Plot 1
    subplot(1,2,1,'ticklabelinterpreter','latex')
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    hold on
    for ii=1:nIter
        nVarLoc=length(optimstruct(ii).population);
        itL=[optimstruct(ii).population(:).objective];
        nVarLoc=length(itL);
        iterRes(ii,1:nVarLoc)=itL;
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    for ii=2:2:nIter
        nVarLoc=length(optimstruct(ii).population);
        lSub1(1)=plot(ones(1,nVarLoc)*(ii+2)/2,iterRes(ii,1:nVarLoc),'b.','markersize',5);
    end
    switch dirOptim
        case 'min'
            minRes=iterRes(1:2:end,1);
        case 'max'
            minRes=iterRes(1:2:end,1);
    end
    
    meanRes=mean(iterRes,2);
    stdRes=std(iterRes,0,2);
    lSub1(2)=plot(1:numel(minRes),minRes,'r-');
    %lSub1(3)=plot(1:2:nIter,meanRes,'color',[0.7 0 0]);
    lSub1(3)=plot([0,nIter],[knownOptim(1) knownOptim(1)],'r--');
    
    legend(lSub1,{'Population',['Root Member'],...
        'Theoretical Optimum'},'Location','NorthEast', 'interpreter','latex');
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('$J(\mathbf{x})$', 'interpreter','latex','fontsize',12)
    
    switch dirOptim
        case 'min'
            testOrder=max(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=min([minRes;knownOptim(1)]);
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
        case 'max'
            testOrder=min(minRes);
            orderSol=ceil(-log10(abs(testOrder)));
            box(3)=floor(testOrder*10^orderSol)*10^(-orderSol);
            testOrder=max([minRes;knownOptim(1)]);
            box(4)=ceil(testOrder*10^orderSol)*10^(-orderSol);
    end
    box(1:2)=[0,ceil(nIter/2)+1];
    
    
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
            errMeasure=-(knownOptim(1)-minRes);
        case 'max'
            errMeasure=knownOptim(1)-minRes;
    end
    
    for ii=1:2:nIter
        nFill(ceil(ii/2))=numel(optimstruct(ii).population(1).fill);
    end
    changeIter=find(nFill(2:end)~=nFill(1:end-1));
    iterList=1:2:nIter;
    changeIter=iterList(changeIter+1);
    changeIter=[1,changeIter,nIter];
    %     for ii=1:nIter
    %         nPop=length(optimstruct(ii).population);
    %         for jj=1:nPop
    %             plot3(ones(size(optimstruct(ii).population(jj).fill))*((ii+1)/2),...
    %                 1:numel(optimstruct(ii).population(jj).fill),...
    %                 optimstruct(ii).population(jj).fill)
    %             hold on
    %         end
    %     end
    for ll=1:(numel(changeIter)-1)
        popNominal=zeros([ceil((changeIter(ll+1)-1-changeIter(ll))/2),numel(optimstruct(changeIter(ll)).population(1).fill)]);
        objNominal=zeros([ceil((changeIter(ll+1)-1-changeIter(ll))/2),1]);
        kk=0;
        for ii=changeIter(ll):2:(changeIter(ll+1)-1)
            kk=kk+1;
            popNominal(kk,:)=optimstruct(ii).population(1).fill;
            objNominal(kk)=optimstruct(ii).population(1).objective;
        end
        
        
        iterVec=((changeIter(ll):2:(changeIter(ll+1)-1))+1)/2;
        varVec=(1:length(popNominal(1,:)))/length(popNominal(1,:));
        
        [iterGrid,varGrid]=meshgrid(iterVec,varVec);
        [objGrid,~]=meshgrid(objNominal,varVec);
        
        surf(iterGrid,varGrid,popNominal',objGrid)
        hold on
    end
    colorbar;
    
    xlabel('Iteration', 'interpreter','latex','fontsize',12)
    ylabel('Variable Index', 'interpreter','latex','fontsize',12)
    zlabel('Variable Value', 'interpreter','latex','fontsize',12)
    set(axh,'ticklabelinterpreter','latex')
    
    
    
    
end

function [h]=PlotAnalyticalPath(paramoptim,optimstruct,objectiveName)
    
    varExtract={'iterGap','knownOptim','desVarRange','optimMethod','nDesVar','direction'};
    [iterGap,knownOptim,desVarRange,optimMethod,nDesVar,dirOptim]=ExtractVariables(varExtract,paramoptim);
    
    h=figure('Name','Optimisation Path','Position',[20 100 1000 600]);
    nIter=numel(optimstruct);
    if CheckIfGradient(optimMethod)
        for ii=1:nIter
            rootPath(ii,1:nDesVar)=(optimstruct(ii).population(1).fill);
            rootDes(ii)=(optimstruct(ii).population(1).objective);
        end
    else
        for ii=1:nIter
            nVarLoc=length(optimstruct(ii).population);
            iterRes(ii,1:nVarLoc)=[optimstruct(ii).population(:).objective];
            lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);
        end
        switch dirOptim
            case 'min'
                [minRes,minPos]=min(iterRes,[],2);
            case 'max'
                [minRes,minPos]=max(iterRes,[],2);
        end
        for ii=1:nIter
            rootPath(ii,1:nDesVar)=optimstruct(ii).population(minPos).fill;
            rootDes(ii)=optimstruct(ii).population(minPos).objective;
        end
    end
    
    for ii=1:2:nDesVar
%         l=plot3(rootPath(:,ii),rootPath(:,ii+1),rootDes,'r-*','linewidth',2);
        l=plot(rootPath(:,ii),rootPath(:,ii+1),'r-*','linewidth',2);
        l.DisplayName=['DesVar : ',int2str([ii ii+1])];
        hold on
    end
    legend
    x=linspace(desVarRange(1),desVarRange(2),200);
    [X1,X2]=meshgrid(x,x);
    Y=zeros(size(X1));
    for ii=1:numel(Y)
        Y(ii)=eval([objectiveName,'Function([X1(ii),X2(ii)])']);
    end
    %surf(X1,X2,Y,log10(Y))
    contour(X1,X2,log10(Y),40)
    
%     for ii=1:numel(optimstruct)
%         pathIter=vertcat(optimstruct(ii).population(:).fill);
%         resIter=[optimstruct(ii).population(:).objective];
%         for jj=1:2:nDesVar
%             %plot3(pathIter(:,jj),pathIter(:,jj+1),resIter,'+')
%             plot(pathIter(:,jj),pathIter(:,jj+1),'+')
%         end
%     end
    
end

function [h]=PlotGradientsRefine(paramoptim,optimstruct)
    
    supportOptim=paramoptim.optim.supportOptim;
    maxNums=numel(supportOptim.hist(end).gradfk);
    for ii=numel(supportOptim.hist)-1:-1:1
        if numel(supportOptim.hist(ii).gradfk)<maxNums
            break
        end
    end
    ii=ii+1;
    iterstart=(ii-1)*2+1;
    supportOptim.hist=supportOptim.hist(ii:end);
    paramoptim.optim.supportOptim=supportOptim;
    if (numel(optimstruct)-iterstart)>=4
        [h,~]=PlotGradients(paramoptim,optimstruct(iterstart:end));
    else
        h=figure;
    end
end

function [h,directionChange]=PlotGradients(paramoptim,optimstruct)
    
    normVec=@(vec) sqrt(sum(vec.^2,2));
    normaliseArray=@(array)  array./(sqrt(sum(array.^2,2))*ones(1,size(array,2))); 
    
    supportOptim=paramoptim.optim.supportOptim;
    
    h=figure('Name',['Gradients_',ExtractVariables({'optimCase'},paramoptim)]...
        ,'Position',[ 100 150 1400 700]);
    symDesVarList=ExtractVariables({'symDesVarList'},paramoptim);
    rmCol=symDesVarList(2,:);
    subplot(2,2,1)
    grads=vertcat(supportOptim.hist(:).gradfk);
    grads(:,rmCol)=[];
    surf(log10(abs(grads)))
    ylabel('iteration')
    xlabel('design variable')
    title('gradients')
    
    view(0,90)
    subplot(2,2,2)
    grads=vertcat(supportOptim.hist(:).prevDir);
    grads(:,rmCol)=[];
    surf(log10(abs(grads)))
    ylabel('iteration')
    xlabel('design variable')
    title('previous direction')
    view(0,90)
    subplot(2,3,4)
    grads=vertcat(supportOptim.hist(:).prevStep);
    surf(((grads)),(abs((grads))))
    ylabel('iteration')
    xlabel('design variable')
    title('previous direction')
    view(0,90)
    subplot(2,3,5)
    hold on
    grads=vertcat(supportOptim.hist(:).gradfk);
    grads(:,rmCol)=[];
    directionChange.gradNorm=normVec(grads);
    gradsm1=vertcat(supportOptim.hist(:).gradfkm1);
    gradsm1(:,rmCol)=[];
    directionChange.Grad=dot(normaliseArray(grads),normaliseArray(gradsm1),2);
    directionChange.Grad(1)=1;
    gradStep=vertcat(supportOptim.hist(:).prevStep);
    gradStep(:,rmCol)=[];
    directionChange.Step=dot(normaliseArray(grads),normaliseArray(gradStep),2);
    gradDir=vertcat(supportOptim.hist(:).prevDir);
    gradDir(:,rmCol)=[];
    directionChange.Dir=dot(normaliseArray(grads),normaliseArray(gradDir),2);
    directionChange.grads=grads;
    directionChange.gradsm1=gradsm1;
    kk=1;
    fillPrec=zeros(size(optimstruct(1).population(1).fill));
    for ii=1:2:numel(optimstruct)
        changePos(kk)=normVec(fillPrec-optimstruct(ii).population(1).fill);
        fillInf(kk,1:numel(optimstruct(ii).population(1).fill))...
            =optimstruct(ii).population(1).fill;
        fillPrec=optimstruct(ii).population(1).fill;
        kk=kk+1;
    end
    changePos(1)=ExtractVariables({'startVol'},paramoptim);
    directionChange.Pos=changePos;
    
    step=1:numel(supportOptim.hist);
    ii=1;
    l(ii)=plot(step,directionChange.Grad);
    l(ii).DisplayName='Change in direction - gradient';
    ii=ii+1;
    l(ii)=plot(step,directionChange.Step);
    l(ii).DisplayName='gradient // direction';
    ii=ii+1;
    l(ii)=plot(step,directionChange.Dir);
    l(ii).DisplayName='gradient // step';
    ii=ii+1;
    try
        scaleEvol=[supportOptim.hist(:).scale];
        l(ii)=plot(step,scaleEvol);
        l(ii).DisplayName='scale';
        ii=ii+1;
    catch
    end
    try
        iterEvol=[supportOptim.hist(:).iter];
        l(ii)=plot(step,iterEvol);
        l(ii).DisplayName='iter since last refresh';
        ii=ii+1;
    catch
    end
    l(ii)=plot(step,changePos(end-numel(step)+1:end));
    l(ii).DisplayName='Length of movement';
    ii=ii+1;
    l(ii)=plot(step,directionChange.gradNorm);
    l(ii).DisplayName='Norm of Gradient';
    ii=ii+1;
    legend(l)
    ax=subplot(2,3,6);
    if isfield(optimstruct(1).population(1).additional,'res')
        [ax]=PlotDiffVsCFDConv(optimstruct,ax);
    end
end

function [ax,h]=PlotDiffVsCFDConv(optimstruct,ax)
    
    if nargin<2
        h=figure;
        ax=axes;
    else
        h=ax.Parent;
    end
    
    hold on
    baseStruct=struct('dat',[],'mean',[],'max',[],'min',[]);
    convTest=struct('objres',repmat(baseStruct,size(1:2:length(optimstruct)))...
        ,'diffres',repmat(baseStruct,size(1:2:length(optimstruct))));
    kk=1;
    convFields=fieldnames(convTest);
    
    for ii=1:2:length(optimstruct)
        keepdiffres=true([1, numel(optimstruct(ii).population)-1]);
        for jj=2:numel(optimstruct(ii).population)
            if ~isempty(optimstruct(ii).population(jj).additional.res)
                convTest.objres(kk).dat(jj-1)=[optimstruct(ii).population(jj).objective]./...
                    10.^[optimstruct(ii).population(jj).additional.res];
                convTest.diffres(kk).dat(jj-1)=abs([optimstruct(ii).population(jj).objective]-...
                    optimstruct(ii).population(1).objective)./...
                    10^[optimstruct(ii).population(jj).additional.res];
            else 
                keepdiffres(jj-1)=false;
            end
        end
        convTest.diffres(kk).dat=convTest.diffres(kk).dat(keepdiffres);
        convTest.objres(kk).dat=convTest.objres(kk).dat(keepdiffres);
        for jj=1:numel(convFields)
            convTest.(convFields{jj})(kk).mean=mean(convTest.(convFields{jj})(kk).dat);
            convTest.(convFields{jj})(kk).max=max(convTest.(convFields{jj})(kk).dat);
            convTest.(convFields{jj})(kk).min=min(convTest.(convFields{jj})(kk).dat);
        end
        kk=kk+1;
    end
    x=1:length(convTest.objres);
    l3=plot(ax,x,[convTest.objres(:).mean],...
       x,[convTest.diffres(:).mean]);
    
    l3(1).DisplayName='obj/res';
    c{1}=l3(1).Color;
    c{2}=l3(2).Color;
    l3(2).DisplayName='diff/res';
    
    l2=plot(ax,x,[convTest.objres(:).min],'--',...
        x,[convTest.diffres(:).min],'--');
    [l2.Color]=deal(c{:});
    l2=plot(ax,x,[convTest.objres(:).max],'--',...
        x,[convTest.diffres(:).max],'--');
    [l2.Color]=deal(c{:});
    legend(l3);
    ax.YScale='log'; 
end


%% Rebuild Test population 

function [teststruct]=RebuildFromIter(iterstruct,out)
    
    teststruct=iterstruct;
    kk=0;
    numel(iterstruct)
    
    for ii=1:numel(out)
        if isdir(out(ii).rootDir)
            [iterDir,iterName]=FindDir(out(ii).rootDir,'iteration_',1);
            for jj=1:numel(iterDir)
                if isdir(iterDir{jj})
                    [iterbin,iterBinName]=FindDir(iterDir{jj},'population_iteration',0);
                    if ~isempty(iterbin)
                        itertest=load(iterbin{1});
                        ll=(regexp(iterBinName{1},'[0-9]*','match'));
                        ll=str2double(ll{1});
                        teststruct(ll).population=itertest.population;
                    else
                        kk=kk+1;
                    end
                else
                    kk=kk+1;
                end
            end
        else
            warning('Non Existing Root Dir')
        end
    end
    
%     for ii=1:numel(iterstruct)
%         curiterpath='';
%         nPop=numel(iterstruct(ii).population);
%         jj=1;
%         while isempty(curiterpath) && jj<=nPop
%             curiterpath=iterstruct(ii).population(jj).location;
%             jj=jj+1;
%         end
%         
%         if ~isempty(curiterpath)
%             curiterpath=regexp(curiterpath,'^.*iteration_[0-9]*','match');
%             if isdir(curiterpath{1})
%             [iterbin,~]=FindDir(curiterpath{1},'population_iteration',0);
%             if ~isempty(iterbin)
%                 itertest=load(iterbin{1});
%                 teststruct(ii).population=itertest.population;
%             else
%                 kk=kk+1;
%             end
%             else
%                 kk=kk+1;
%             end
%             
%         else
%             kk=kk+1;
%         end
%     end 
    
    if kk>0
        warning([int2str(kk),' test populations not found, live population will be displayed in their place'])
    end
    
end

%% Binaries

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

function []=GenerateIterResultBinary(resultDirectory,marker,optimstruct,paramoptim)
    
    fileName=[resultDirectory,'\OptimRes_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'optimstruct');
    fileName=[resultDirectory,'\FinalParam_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'paramoptim');
    
end

function []=GenerateTestResultBinary(resultDirectory,marker,optimstruct,paramoptim)
    
    fileName=[resultDirectory,'\OptimTestRes_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'optimstruct');
    
end

function [dat]=GenerateOptimalSolDir(resultDirectory,markerSmall,optimDirection,optimstruct)
    
    [~,posOpt]=eval([optimDirection,'([optimstruct(end).population(:).objective])']);
    optimsolution=optimstruct(end).population(posOpt);
    resultDirectory=[resultDirectory,'\Optimal_',markerSmall];
    resultDirectory=MakePathCompliant(resultDirectory);
   
    profileDir=optimsolution.location;
    
    %copyfile(profileDir,resultDirectory,'f');
    compType=computer;
    if strcmp(compType(1:2),'PC')
        copyfileRobust(profileDir,resultDirectory);
    else
        [~,~]=system(['rsync -r --exclude=*.*[0-9]* ''',profileDir,'/'' ''',resultDirectory,'''']);
        
    end
    c=dir(resultDirectory);
    isFileName=false;
    ii=0;
    while(~isFileName) && ii<numel(c)
        ii=ii+1;
        isFileName=~isempty(regexp(c(ii).name,'restart', 'once'));
    end
    if isFileName
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
    else
        dat.A=0;
        coord=[0,0];
        dat.xMin=min(coord(:,1));
        dat.xMax=max(coord(:,1));
        dat.t=max(coord(:,2))-min(coord(:,2));
        dat.nPoints=length(coord(:,1));
    end
    
end


