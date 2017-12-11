function []=TecOutErrorWithProfile()
    
    
end

%% Create the additional tecsubfile with error as vel

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

%% 