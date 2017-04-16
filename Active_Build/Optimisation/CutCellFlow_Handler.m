

function [obj]=CutCellFlow_Handler(paramoptim,boundaryLoc)
    
    % copy standard executables and setting
    
    
    [targFolder]=PrepareCFDFolder(paramoptim,boundaryLoc);
    GenerateMesh(paramoptim,targFolder);
    [obj]=SolveFlow(paramoptim,targFolder);
end


%% Main Handler blocks

function [targFolder]=PrepareCFDFolder(paramoptim,boundaryLoc)
    varExtract={'CFDfolder','nMach'};
    [CFDfolder,nMach]=ExtractVariables(varExtract,paramoptim);
    
    compType=computer;
    boundaryLoc=MakePathCompliant(boundaryLoc);
    CFDfolder=MakePathCompliant(CFDfolder);
    
    inDir=what(boundaryLoc);
    boundaryLoc=inDir.path;
    targFolder=[boundaryLoc,filesep,'CFD'];
    targFolder=MakePathCompliant(targFolder);
    if strcmp(compType(1:2),'PC')
        copyfile(CFDfolder,targFolder)
    else
        flowCommand=['cp -rp ''',CFDfolder,''' ''',targFolder,''''];
        [status,stdout]=system(flowCommand);
    end
    if numel(nMach)<3
        nMach=[1 nMach(1) 0.25];
    end
    RestartModifiedSettings(targFolder,'mach',1,nMach);
    CopyBoundaryFile(boundaryLoc,targFolder);
    
end


function []=GenerateMesh(paramoptim,targFolder)
    varExtract={'isSymFlow','parentMesh'};
    [isSymFlow,parentMesh]=ExtractVariables(varExtract,paramoptim);
    compType=computer;
    if strcmp(compType(1:2),'PC')
        extStr='bat';
    else
        extStr='sh';
    end
    
    % Either Generate or Load Mesh
    parentMesh=[parentMesh,filesep,'CFD',filesep,'griduns'];
    isLoad=numel(dir(parentMesh))>0;
    if isLoad
        meshGenCommand=['cp "',parentMesh,'" "',targFolder,filesep,'griduns"'];
        isDeformation=true;
    else
        meshGenCommand=['"',targFolder,filesep,'RunCartCell.',extStr,'"'];
        isDeformation=false; % probably the right choice
    end
    [status,stdout]=system(meshGenCommand);
    
    % Mesh Symmetry
    if isSymFlow
        meshSymCommand=['"',targFolder,filesep,'CartCellSym.',extStr,'"'];
        [status,stdout]=system(meshSymCommand);
    end
    % Mesh Deformation
    if isDeformation
        meshDeformCommand=['"',targFolder,filesep,'gridtoxyz.',extStr,'"'];
        [status,stdout]=system(meshDeformCommand);
    end
    % Check for errors in deformation and revert to remeshing
    
    
    
    
end


function [obj]=SolveFlow(paramoptim,targFolder)
    
    varExtract={'stoponerror','targConv','restartIter','lengthConvTest',...
        'maxRestart','startIterFlow','flowRestart','parentMesh'};
    [stoponerror,targConv,restartIter,lengthConvTest,maxRestart,...
        startIterFlow,flowRestart,parentMesh]...
        =ExtractVariables(varExtract,paramoptim);
    
    
    if flowRestart && (numel(parentMesh)>0)
        flowPath=FindDir([parentMesh,filesep,'CFD'],'flow.dump',false);
        if numel(flowPath)>0
            flowCopy=['cp "',flowPath{1},'" "',targFolder,filesep,'flow.dump"'];
            [status,stdout]=system(flowCopy);
            RestartModifiedSettings(targFolder,'restart',12,1);
        else
            flowRestart=false;
        end
    end
    
    RestartModifiedSettings(targFolder,'iter',7,...
        [2, startIterFlow, targConv]);
    compType=computer;
    errFlag=true;%CutCellErrorHandling(endstr,false);
    kk=0;
    iterN=0;
    while sum(errFlag) && kk<10
        endstr=RunFlowSolverOnly(compType,targFolder);
        errFlag=CutCellErrorHandling(endstr,false);
        RestartModifiedSettings(targFolder,'cfl',8)
        kk=kk+1;
        if kk==1 && flowRestart && sum(errFlag)
            RestartModifiedSettings(targFolder,'restart',12,0);
            warning('Flow restart failed with the following message: %s',endstr)
        end
    end
    
    if kk<10
        [obj]=ExtractFinalData(targFolder,iterN,sum(errFlag));
        iterN=obj.iter;
        % test convergence
        
        [isfinished,restartNormal,restartCFL,theoretConvIter]=...
            TestCutCellConvergence(obj.res,targConv,targFolder,lengthConvTest);
        kk=0;
        while ~isfinished && kk<maxRestart
            theoretConvIter(theoretConvIter<=0)=Inf;
            theoretConvIter(theoretConvIter<=20)=20;
            if restartNormal || restartCFL
                RestartModifiedSettings(targFolder,'restart',12,1);
                RestartModifiedSettings(targFolder,'iter',7,...
                    [2, min([restartIter,theoretConvIter]), targConv]);
            end
            
            if restartCFL
                RestartModifiedSettings(targFolder,'cfl',8);
            end
            
            endstr=RunFlowSolverOnly(compType,targFolder);
            kk=kk+1;
            
            if sum(errFlag)
                isfinished=false;
                restartNormal=false;
                restartCFL=true;
            else
                [obj]=ExtractFinalData(targFolder,iterN,sum(errFlag));
                [isfinished,restartNormal,restartCFL,theoretConvIter]=...
                    TestCutCellConvergence(obj.res,targConv,targFolder,lengthConvTest);
                iterN=obj.iter;
            end
        end
        errFlag=CutCellErrorHandling(endstr,stoponerror);
    end
    
    [obj]=ExtractFinalData(targFolder,iterN,sum(errFlag));
    
end

%% Flow solver operations

function [obj]=ExtractFinalData(targFolder,iter,isErr)
    
    resFile=[targFolder,filesep,'res_hist.dat'];
    
    fid=fopen(resFile,'r');
    structName={'iter','res','cl','cm','cd','cx','cy'};
    
    if ~isErr
        fseek(fid,-300,1); % Seek the end of file
        while ~feof(fid)
            str=fgetl(fid);
        end
        data=str2num(str); %#ok<ST2NM>
        for ii=1:length(structName)
            obj.(structName{ii})=data(ii);
        end
    else
        for ii=1:length(structName)
            obj.(structName{ii})=1000;
        end
    end
    fclose('all');
    obj.iter=obj.iter+iter;
    
end

function endstr=RunFlowSolverOnly(compType,targFolder)
    
    if strcmp(compType(1:2),'PC')
        extStr='bat';
    else
        extStr='sh';
    end
    
    flowCommand=['"',targFolder,filesep,'RunEulerFlow.',extStr,'"'];
    [status,stdout]=system(flowCommand);
    endstr=stdout(end-200:end);
    
end

function []=CopyBoundaryFile(boundaryLoc,targFolder)
    
    c=dir(boundaryLoc);
    kk=1;
    nDir=length(c);
    flag=false;
    while ~flag && kk<=nDir
        flag=~isempty(regexp(c(kk).name,'boundary', 'once'));
        kk=kk+1;
    end
    
    if kk==nDir+1 && ~flag
        error('Boundary file not saved, flow cannot be computed')
    end
    origFile=[boundaryLoc,filesep,c(kk-1).name];
    targFile=[targFolder,filesep,'boundary.dat'];
    copyfile(origFile,targFile);
    dispPath=FindDir(boundaryLoc,'displacements',false);
    if numel(dispPath)>0
        targDisp=[targFolder,filesep,'displacements.xyz'];
        copyfile(dispPath{1},targDisp);
        dispPath=FindDir(boundaryLoc,'surface',false);
        targDisp=[targFolder,filesep,'surface.xyz'];
        copyfile(dispPath{1},targDisp);
    end
end

%% Error handling

function [errFlag]=CutCellErrorHandling(endstr,stoponerror)
    
    errorTerms={'exception','neg','griduns','error','IEEE_INVALID_FLAG','IEEE_UNDERFLOW_FLAG','lam'};
    [errFlag,errorstr]=CutCellErrorDetection(errorTerms,endstr);
    
    if sum(errFlag) && stoponerror
        error(errorstr)
    end
    
end

function [errFlag,errorstr]=CutCellErrorDetection(errorTerms,endstr)
    
    errorList=regexp(endstr,errorTerms, 'once');
    errorstr='';
    errFlag=false([1 length(errorTerms)]);
    for ii=1:length(errorTerms)
        if ~isempty(errorList{ii})
            errorstr=char(errorstr,['Error detected: ',errorTerms{ii}]);
            errFlag(ii)=true;
        end
    end
    
end

%% Convergence Issue Handling

function [isfinished,restartNormal,restartCFL,theoretConvIter]=...
        TestCutCellConvergence(finConv,targConv,targFolder,numAverage)
    
    isfinished=(finConv<targConv) || (finConv==0);
    theoretConvIter=0;
    restartNormal=false;
    restartCFL=false;
    if ~isfinished
        [resArray]=ExtractConvergenceRate(targFolder,numAverage);
        [isConverging,theoretConvIter]=TestConvergenceRate(resArray,targConv);
        if isConverging
            restartNormal=true;
        else
            restartCFL=true;
        end
    end
    
    
end

function [streamArray]=ExtractConvergenceRate(targFolder,numAverage)
    
    fidHist=fopen([targFolder,filesep,'res_hist.dat'],'r');
    
    % skip header
    hLines=2;
    for ii=1:hLines
        fgetl(fidHist);
    end
    
    streamArray=zeros([numAverage,2]);
    kk=0;
    while(~feof(fidHist))
        ii=mod(kk,numAverage)+1;
        actLine=str2num(fgetl(fidHist)); %#ok<ST2NM>
        streamArray(ii,:)=actLine(1:2);
        kk=kk+1;
    end
    
    reorderStream=[ii+1:min([kk,numAverage]),1:ii];
    streamArray=streamArray(reorderStream,:);
    fclose('all');
    
    
end

function [isConverging,theoretConvIter]=TestConvergenceRate(resArray,targConv)
    isConverging=true;
    meanRes=mean(resArray(:,2));
    stdRes=std(resArray(:,2));
    minRes=min(resArray(:,2));
    finRes=resArray(end,2);
    startRes=resArray(end,2);
    lastIter=max(resArray(:,1));
    
    matLin=[resArray(:,1),ones(size(resArray(:,1)))];
    coeffLin=(matLin'*matLin)\matLin'*resArray(:,2);
    linRMS=sqrt(mean((matLin*coeffLin).^2));
    theoretConvIter=round((targConv-coeffLin(2))/coeffLin(1)-lastIter);
    
    if theoretConvIter<0
        isConverging=false;
    end
    
    
end

%% Settings file rewrite

function []=RestartModifiedSettings(targFolder,typeChange,lineNum,varargin)
    
    origSettings=[targFolder,filesep,'settings'];
    saveSettings=[targFolder,filesep,'settings_archive'];
    
    copyfile(origSettings,saveSettings);
    
    fidSetW=fopen(origSettings,'w');
    fidSetR=fopen(saveSettings,'r');
    
    for ii=1:lineNum-1
        fprintf(fidSetW,fgets(fidSetR));
    end
    
    cflLine=fgetl(fidSetR);
    switch typeChange
        case 'cfl'
            cfloutstr=ReplaceCFLNum(cflLine);
        case 'iter'
            cfloutstr=ReplaceIter(cflLine,varargin{1});
        case 'restart'
            cfloutstr=ReplaceRestart(cflLine,varargin{1});
        case 'mach'
            cfloutstr=ReplaceMachNum(cflLine,varargin{1});
    end
    
    fprintf(fidSetW,cfloutstr);
    
    while (~feof(fidSetR))
        fprintf(fidSetW,fgets(fidSetR));
    end
    fclose('all');
end

function strOut=ReplaceCFLNum(strIn)
    
    cflStr=regexp(strIn,'\d*\.\d*','match','once');
    strOut=regexprep(strIn,cflStr,'%f');
    cflNum=str2double(cflStr);
    cflNum=cflNum/2;
    strOut=sprintf([strOut,'\n'],cflNum);
    
end

function strOut=ReplaceMachNum(strIn,nMach)
    
    cflStr=regexp(strIn,'\d*\s*\d*\.\d*\s*\d*\.\d*','match','once');
    strOut=regexprep(strIn,cflStr,'%i %.2f %.2f');
    
    strOut=sprintf([strOut,'\n'],nMach(1),nMach(2),nMach(3));
    
end

function strOut=ReplaceIter(strIn,iterNum)
    
    cflStr=regexp(strIn,'\d*\s*\d*\s*\d*\.\d*','match','once');
    if isempty(cflStr)
        cflStr=regexp(strIn,'\d*\s*\d*\s*-\d*\.\d*','match','once');
    end
    if isempty(cflStr)
        cflStr=regexp(strIn,'\d*\s*\d*\s*\d*','match','once');
    end
    if isempty(cflStr)
        cflStr=regexp(strIn,'\d*\s*\d*\s*-\d*','match','once');
    end
    strOut=regexprep(strIn,cflStr,'%i %i %.2f ');
    
    strOut=sprintf(['%i %i %.2f ','\n'],iterNum(1),iterNum(2),iterNum(3));
    
end

function strOut=ReplaceRestart(strIn,isrestart)
    
    cflStr=regexp(strIn,'\d*');
    strOut=strIn;
    strOut=[strOut(1:cflStr(2)),' ',strOut((cflStr(2)+1):end)];
    strOut(cflStr(2):cflStr(2)+1)='%i';
    strOut=sprintf([strOut,'\n'],isrestart);
    
end
