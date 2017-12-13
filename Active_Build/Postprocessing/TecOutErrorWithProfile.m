function []=TecOutErrorWithProfile(pathStr)
    
    [pathstruct,optimstruct]=FindRestartsTecsubfile(pathStr);
    for ii=1:numel(pathstruct)
        %[paramoptim]=BuildErrorTecSubfile(pathstruct(ii));
    end
    load(pathstruct(end).param)
    marker='RefinementSummary';
    
    tecPlotFile{1}=['Tec360plt_Flow_',marker,'.plt'];
    tecPlotFile{2}=['Tec360plt_Snak_',marker,'.plt'];
    tecPlotPre{1}=['Tec360plt_Flow_',marker,'_pre.plt'];
    tecPlotPre{2}=['Tec360plt_Snak_',marker,'_pre.plt'];
    
    [FID]=OpenOptimumSnakLayFile(pathStr,marker);
    PersnaliseLayFile(FID,tecPlotPre(2:-1:1));
    tecPlotFile{1}=[pathStr,filesep,tecPlotFile{1}];
    tecPlotFile{2}=[pathStr,filesep,tecPlotFile{2}];
    axisRatio=ExtractVariables({'axisRatio'},paramoptim.parametrisation);
    ExtractOptimalSnake(optimstruct,pathStr,'min',...
            tecPlotFile,axisRatio,paramoptim,{pathstruct(2:end).dir});
end

%% Create the additional tecsubfile with error as vel

function [pathstruct,optimstruct]=FindRestartsTecsubfile(pathStr)
    
    dirPath=FindDir(pathStr,'Dir_',1);
    
    % Reorder paths to be in order of refinement
    dirOrd=regexp(dirPath,'[0-9]*$','match');
    dirOrd=[dirOrd{:}];
    dirPath(cellfun(@isempty,dirOrd))=[];
    dirOrd(cellfun(@isempty,dirOrd))=[];
    dirOrd=cellfun(@str2num,dirOrd);
    [~,dirOrd]=sort(dirOrd);
    dirPath=dirPath(dirOrd);
    
    for ii=1:numel(dirPath)+1
        if ii==1
            optimPath={MakePathCompliant([dirPath{max(ii-1,1)},...
                '\iteration_1\'])};
        else
            optimPath=FindDir(dirPath{max(ii-1,1)},'Optimal_',1);
        end
        restartTemp=FindDir(optimPath{1},'restart_',0);
        [tecTemp,tecName]=FindDir(optimPath{1},'tecsubfile',0);
        tt=find(~cellfun(@isempty,regexp(tecName,'^tecsubfile_.*\.dat$')));
        pathstruct(ii).restart=restartTemp{1};
        pathstruct(ii).tecplot=tecTemp{tt};
        pathstruct(ii).param=FindDir(dirPath{max(ii-1,1)},'FinalParam',0);
        pathstruct(ii).param=pathstruct(ii).param{1};
        pathstruct(ii).dir=dirPath{max(ii-1,1)};
        
        optimstruct(ii).population=struct('objective',1,'location',optimPath{1});
    end
    
end

function [paramoptim]=BuildErrorTecSubfile(pathStruct)
    
    load(pathStruct.param)
    load(pathStruct.restart)
    paramoptim=SetVariables({'profileComp'},{'normdist'},paramoptim);
    [errorMeasure,h,targCoord,analysisLoop]=InverseDesign_ErrorTopo(paramoptim,loop);
    close(h)
    snaxel=repmat(struct('index',0,'snaxnext',0,...
    'v',0),[sum([analysisLoop(:).nPts]) 1]);
    snakposition=repmat(struct('index',0,'coord',[0 0]...
        ,'vector',[0 0]),[sum([analysisLoop(:).nPts]) 1]);
    
    kk=1;
    for ii=1:numel(analysisLoop)
        
        kkStart=kk;
        for jj=1:analysisLoop(ii).nPts
            snaxel(kk).index=kk;
            snaxel(kk).snaxnext=kk+1;
            snaxel(kk).v=analysisLoop(ii).localerror(jj);
            snakposition(kk).index=kk;
            snakposition(kk).coord=analysisLoop(ii).coord(jj,:);
            kk=kk+1;
        end
        snaxel(kk-1).snaxnext=kkStart;
    end
    tCell=regexp(regexp(pathStruct.restart,'restart_.*$','match'),'[0-9]*','match');
    time=str2num([tCell{1}{1},'.',tCell{1}{2}]);
    [cellMesh]=SnaxelToCellOut(snaxel,snakposition,10,time);
    FID=fopen(pathStruct.tecplot,'a');
    WriteToFile(cellMesh,FID)
    fclose(FID);
end

% Function to concatenate it all.

%% Optimisation Output

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
        val={ii,[3 4 10],minIterRootDirNum(ii)*5,...
            minIterRootDirNum(ii)*5};
        nOccur=[3 3 1 1];
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

function [cellMesh]=SnaxelToCellOut(snaxel,snakposition,strandID,time)
    
    
    coordDat=vertcat(snakposition(:).coord);
    vectorDat=vertcat(snakposition(:).vector);
    velDat=[snaxel(:).v]';
    %velDat=[[snaxel(:).isfreeze]']*2+1;
    %velDat=[[snaxel(:).orderedge]'];
    vectorDat(:,1)=vectorDat(:,1).*velDat;
    vectorDat(:,2)=vectorDat(:,2).*velDat;
    vectorDat=[vectorDat,velDat];
    vertIndex=[snakposition(:).index];
    connDat=[[snaxel(:).index]',[snaxel(:).snaxnext]'];
    
    [cellMesh]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID,time);
    
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


%% 