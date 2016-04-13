

function [obj]=CutCellFlow_Handler(paramoptim,boundaryLoc)
    
    % copy standard executables and setting
    varExtract={'CFDfolder','stoponerror','targConv','restartIter','lengthConvTest',...
        'maxRestart'};
    [CFDfolder,stoponerror,targConv,restartIter,lengthConvTest,maxRestart]...
        =ExtractVariables(varExtract,paramoptim);
    inDir=what(boundaryLoc);
    boundaryLoc=inDir.path;
    targFolder=[boundaryLoc,filesep,'CFD'];
    copyfile(CFDfolder,targFolder)
    CopyBoundaryFile(boundaryLoc,targFolder);
    % Run system command
    compType=computer;
    if strcmp(compType(1:2),'PC')
        flowCommand=['"',targFolder,filesep,'RunFlow.bat"'];
        [status,stdout]=system(flowCommand);
    else
        flowCommand=['''',targFolder,filesep,'RunFlow.sh'''];
        [status,stdout]=system(flowCommand);
    end
    
    endstr=stdout(end-200:end);
    errFlag=CutCellErrorHandling(endstr,stoponerror);
    
    while sum(errFlag)
        RestartModifiedSettings(targFolder,'cfl',8)
        endstr=RunFlowSolverOnly(compType,targFolder);
        errFlag=CutCellErrorHandling(endstr,stoponerror);
    end
    
    datastr=endstr(end-60:end);
    data=str2num(datastr); %#ok<ST2NM>
    obj.iter=data(1);
    obj.err=data(2);
    obj.cl=data(3);
    obj.cd=data(4);
    iterN=obj.iter;
    % test convergence
    
    [isfinished,restartNormal,restartCFL,theoretConvIter]=...
        TestCutCellConvergence(obj.err,targConv,targFolder,lengthConvTest);
    kk=0;
    while ~isfinished && kk<maxRestart
        theoretConvIter(theoretConvIter<=0)=Inf;
        theoretConvIter(theoretConvIter<=20)=20;
        if restartNormal || restartCFL
            RestartModifiedSettings(targFolder,'restart',12,1);
            RestartModifiedSettings(targFolder,'iter',7,...
                min([restartIter,theoretConvIter]));
        end
        
        if restartCFL
            RestartModifiedSettings(targFolder,'cfl',8);
        end
        
        endstr=RunFlowSolverOnly(compType,targFolder);
        kk=kk+1;
        errFlag=CutCellErrorHandling(endstr,stoponerror);
        if sum(errFlag)
            isfinished=false;
            restartNormal=false;
            restartCFL=true;
        else
            datastr=endstr(end-60:end);
            data=str2num(datastr); %#ok<ST2NM>
            obj.iter=data(1);
            obj.err=data(2);
            obj.cl=data(3);
            obj.cd=data(4);
            [isfinished,restartNormal,restartCFL,theoretConvIter]=...
                TestCutCellConvergence(obj.err,targConv,targFolder,lengthConvTest);
            iterN=iterN+obj.iter;
        end
    end
    [obj]=ExtractFinalData(targFolder,iterN);
end

function [obj]=ExtractFinalData(targFolder,iter)
    
    resFile=[targFolder,filesep,'res_hist.dat'];
    
    fid=fopen(resFile,'r');
    fseek(fid,-200,1); % Seek the end of file
    fgetl(fid);
    data=num2str(fgetl(fid));
    
    structName={'iter','res','cl','cm','cd','cx','cy'};
    for ii=1:length(structName)
        obj.(structName{ii})=data(ii);
    end
    
    obj.iter=iter;
    
end

function endstr=RunFlowSolverOnly(compType,targFolder)
    
    if strcmp(compType(1:2),'PC')
        flowCommand=['"',targFolder,filesep,'RunOnlyFlow.bat"'];
        [status,stdout]=system(flowCommand);
    else
        flowCommand=['''',targFolder,filesep,'RunOnlyFlow.sh'''];
        [status,stdout]=system(flowCommand);
    end
    
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
    
end

%% Error handling

function [errFlag]=CutCellErrorHandling(endstr,stoponerror)
    
    errorTerms={'exception','neg','griduns','error'};
    [errFlag]=CutCellErrorDetection(errorTerms,endstr);
    
    if sum(errFlag) && stoponerror
        error(errorstr)
    end
    
end

function [errFlag]=CutCellErrorDetection(errorTerms,endstr)
    
    errorList=regexp(endstr,errorTerms, 'once');
    errorstr=[];
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

function strOut=ReplaceIter(strIn,iterNum)
    
    cflStr=regexp(strIn,'\d*','match','once');
    strOut=regexprep(strIn,cflStr,'%i');
    
    strOut=sprintf([strOut,'\n'],iterNum);
    
end

function strOut=ReplaceRestart(strIn,isrestart)
    
    cflStr=regexp(strIn,'\d*');
    strOut=strIn;
    strOut=[strOut(1:cflStr(2)),' ',strOut((cflStr(2)+1):end)];
    strOut(cflStr(2):cflStr(2)+1)='%i';
    strOut=sprintf([strOut,'\n'],isrestart);
    
end