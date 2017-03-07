function [] = include_PostProcessing()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%% General
function []=savefig(figh,figname)

print(figh,figname,'-djpeg','-r600')

end

function [marker,t]=GenerateResultMarker(typDat)
    
    t=now;
    marker=[datestr(t,'yyyy-mm-ddTHHMMSS')...
        ,'_',typDat];
    
end

function [resultDirectory]=GenerateResultDirectoryName(marker,resultRoot,...
        archiveName,t)
    if ~exist('t','var'),t=now;end
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,filesep,archiveName,filesep,dateSubFolders,...
        filesep,'Dir_',marker];
    
    resultDirectory=MakePathCompliant(resultDirectory);
    
    mkdir(resultDirectory)
end

function pathName=MakePathCompliant(pathName)
    
    compStr=computer;
    if strcmp(compStr(1:2),'PC')
        pathName=regexprep(pathName,'/','\\');
    else
        
        pathName=regexprep(pathName,'\\','/');
    end
end

function []=WriteToFile(cellLoops,FID)
    % writes the data to the file
    
    
    for ii=1:length(cellLoops)
        if numel(cellLoops{ii})>0
            for jj=1:length(cellLoops{ii}(:,1))
                fprintf(FID,'%s \n',cellLoops{ii}(jj,:));
            end
        else
            fprintf(FID,'\n');
        end
    end
end

function []=ReorganiseSubPlots(ax,sizSubPlot,outerPad,interPad,fontSizes)
    
    if ~exist('outerPad','var'),outerPad=[0.1,0.1,0.05,0.05];end
    if ~exist('interPad','var'),interPad=[0.05,0.05];end
    if ~exist('fontSizes','var'),fontSizes=[10,12];end
    
    xHeight=1-(outerPad(1)+outerPad(2));
    xLength=(xHeight-interPad(1)*(sizSubPlot(1)-1))/sizSubPlot(1);
    xStart=outerPad(1);
    
    yHeight=1-(outerPad(3)+outerPad(4));
    yLength=(yHeight-interPad(2)*(sizSubPlot(2)-1))/sizSubPlot(2);
    yStart=1-outerPad(4)-yLength;
    
    for ii=1:length(ax)
        jj=ii;
        xInd=rem(jj-1,sizSubPlot(1));
        yInd=ceil(jj/sizSubPlot(1))-1;
        
        loXCorn=(xInd)*(xLength+interPad(1))+xStart;
        loYCorn=-(yInd)*(yLength+interPad(2))+yStart;
        
        posVec=[loXCorn,loYCorn,xLength,yLength];
        set(ax(ii),'position',posVec)
        set(ax(ii),'fontsize',fontSizes(1))
        set(get(ax(ii),'xlabel'),'fontsize',fontSizes(2))
        set(get(ax(ii),'ylabel'),'fontsize',fontSizes(2))
        set(get(ax(ii),'zlabel'),'fontsize',fontSizes(2))
        
    end
    
end

function []=CreateValidFolder(pathName)
    error('replace this shit by mkdir')
    
%     if strcmp(compStr(1:2),'PC')
%         pathName=regexprep(pathName,'/','\\');
%         system(['md "',pathName,'"']);
%     else
%         
%         pathName=regexprep(pathName,'\\','/');
%         system(['mkdir ''',pathName,''''])
%     end
end

%% Plotting 

function [datCol]=ProjectColormap(cMap,cDat,cBounds)
    
    if nargin==2
        cBounds=[min(cDat),max(cDat)];
    end
    if cBounds(1)==cBounds(2)
        cBounds(2)=1+cBounds(2);
    end
    nCol=length(cMap(:,1));
    
    cDat(cDat>max(cBounds))=max(cBounds);
    cDat(cDat<min(cBounds))=min(cBounds);
    
    datMap=linspace(cBounds(1),cBounds(2),nCol)';
    
    
    datCol=interp1(datMap,cMap,cDat);
    
end

%% Parameter File

function []=GenerateParameterFile(FID,param,t,marker)
    
    paramCell{1}='% Parameter File';
    paramCell{2}=['% ',datestr(t)];
    paramCell{3}=['% ',marker];
    for ii=1:length(param.structdat.vars)
        paramCell=[paramCell,ExtractVariablePathAndValue2(param,ii)];
        
    end
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function [paramStr]=ExtractVariablePathAndValue(param,varNum)
    
    varstruct=param.structdat.vars(varNum);
    actstruct=param;
    pathVar='param';
    for jj=1:length(varstruct.vec)
        pathVar=[pathVar,'.',param.structdat.fields{varstruct.vec(jj)}];
        actstruct=actstruct.(param.structdat.fields{varstruct.vec(jj)});
    end
    
    varVal=actstruct;
    [varStr]=GenerateVariableString(varVal);
    paramStr=[pathVar,' = ',varStr,';'];
end

function [paramStr]=ExtractVariablePathAndValue2(param,varNum)
    
    varstruct=param.structdat.vars(varNum);
    actstruct=param;
    pathVar{1}='param';
    
    
    for jj=1:length(varstruct.vec)
        tempPath=[];
        tempPath{numel(pathVar)}=[];
        for kk=1:numel(pathVar)
            pathVar{kk}=[pathVar{kk},'.',param.structdat.fields{varstruct.vec(jj)}];
            nTestVar=numel(eval(pathVar{kk}));
            if (jj<length(varstruct.vec)) && nTestVar>1
                for ii=1:nTestVar
                    tempPath{kk}{ii}=[pathVar{kk},'(',int2str(ii),')'];
                end
            else
                tempPath{kk}{1}=pathVar{kk};
            end
        end
        pathVar=[tempPath{:}];
    end
    paramStr{numel(pathVar)}=[];
    for ii=1:numel(pathVar)
        [varStr]=GenerateVariableString(eval(pathVar{ii}));
        paramStr{ii}=[pathVar{ii},' = ',varStr,';'];
    end
end

function [varStr]=GenerateVariableString(startVar)
    
    classVar=class(startVar);
    [m,n]=size(startVar);
    varStr='';
    switch classVar
        case 'char'
            
            openStr='[';
            closeStr=']';
            if ~isempty(startVar)
                for ii=1:m
                    varStrCell{ii,1}=['''',startVar(ii,:),''''];
                end
                [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,1);
                else
                varStr='''''';
            end
        case 'cell'
            
            openStr='{';
            closeStr='}';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=GenerateVariableString(startVar{ii,jj});
                end
            end
            if exist('varStrCell','var')
                [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
            else
                [varStr]=[openStr,closeStr];
            end
        case 'double'
            
            openStr='[';
            closeStr=']';
            [varStrCell{1:max([m 1]),1:max([n 1])}]=deal(' ');
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=num2str(startVar(ii,jj),24);
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        case 'logical'
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    if startVar(ii,jj) 
                        curStr='true';
                    else
                        curStr='false';
                    end
                    varStrCell{ii,jj}=curStr;
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        otherwise
            if ~isempty(regexp(classVar,'int','once'))
                
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=int2str(startVar(ii,jj));
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
            end
            warning('Class is not catered for and will not be printed correctly to parameter file')
    end
    
end

function [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n)

    
    varStr=openStr;
    if m==0
        m=1;
    end
    if n==0
        n=1;
    end
    for ii=1:m-1
        for jj=1:n-1
            varStr=[varStr,varStrCell{ii,jj},','];
        end
        varStr=[varStr,varStrCell{ii,n},';'];
    end
    for jj=1:n-1
        varStr=[varStr,varStrCell{m,jj},','];
    end
    
    varStr=[varStr,varStrCell{m,n},closeStr];
end

function []=CopyDiary(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['diary_',marker,'.log'];
    originalLayFile=[cd,'\Result_Template\Latest_Diary.log'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfile(originalLayFile,MakePathCompliant([writeDirectory,filesep,fileName]))
    
end



%% Boundary Output to .dat file

function []=BoundaryOutput(loop,FID,typeLoop,buildInternal)
    
    if nargin<4
        buildInternal=false;
        if nargin<3
            typeLoop='subdivision';
        end
    end
    if numel(loop)==0
        error('Invalid loop input, check ''fill''')
    end
    % trim loops and extract data
    if ~buildInternal
        [isInternal]=FindInternalLoop(loop);
        loop=loop(~isInternal);
    end
    if numel(loop)==0
        error('Loop is Empty after internal Trimming, check ''fill''')
    end
    loopout=TrimLoops(loop,typeLoop);
    % format numeric data to printable strings
    cellLoops=DataToString(loopout);
    
    % print string data to file
    WriteToFile(cellLoops,FID)
    
end

function [G]=BuildMatlabGeometryMatrix(loop,typeLoop)
    [isInternal]=FindInternalLoop(loop);
    
    nPts=0;
    for ii=1:numel(loop)
        nPts=nPts+size(loop(ii).(typeLoop),1);
    end
    
    G=zeros(7,nPts);
    nDone=0;
    for ii=1:numel(loop)
        nDoneNew=nDone+size(loop(ii).(typeLoop),1);
        padOnes=ones([1,size(loop(ii).(typeLoop),1)]);
        actCoord=loop(ii).(typeLoop);
        if ~loop(ii).isccw
            actCoord=flip(actCoord);
        end
        actCoord=actCoord';
        G(:,nDone+1:nDoneNew)=[padOnes*2;
            actCoord(1,:);
            actCoord(1,[2:end,1]);
            actCoord(2,:);
            actCoord(2,[2:end,1]);
            padOnes*mod(isInternal(ii)+1,2);
            padOnes*mod(isInternal(ii),2)];
        nDone=nDoneNew;
    end
    
end


function []=OutputFreefempp(loop,typeLoop)
    [isInternal]=FindInternalLoop(loop);
    
    
    
end

function [loopout]=TrimLoops(loop,typeLoop)
    % function extracting the data that must be written to the boundary.dat
    % file
    if nargin<2
        typeLoop='subdivision';
    end
    
    nLoop=length(loop);
    
    for ii=1:nLoop
        endInd=length(loop(ii).(typeLoop)(:,1));
        
        isNeg2=sum(sum(loop(ii).(typeLoop)(1:2,:)==loop(ii).(typeLoop)(end-1:end,:)))...
            ==numel(loop(ii).(typeLoop)(1:2,:));
        isNeg1=sum(loop(ii).(typeLoop)(1,:)==loop(ii).(typeLoop)(end,:))...
            ==numel(loop(ii).(typeLoop)(1,:));
        
        if isNeg2
            endInd=endInd-2;
        elseif isNeg1
            endInd=endInd-1;
        else
        
        end
        
        if loop(ii).isccw
            loopout.surf(ii).coord=loop(ii).(typeLoop)(1:endInd,:);
        else
            loopout.surf(ii).coord=loop(ii).(typeLoop)(endInd:-1:1,:);
        end
        loopout.surf(ii).nvertex=size(loopout.surf(ii).coord,1);
        loopout.surf(ii).nfaces=loopout.surf(ii).nvertex;
    end
    
    loopout.total.nvertex=sum([loopout.surf(:).nvertex]);
    loopout.total.nfaces=sum([loopout.surf(:).nfaces]);
    loopout.total.nloops=nLoop;
    
end

function cellLoops=DataToString(loopout)
    % Transforms the data organised in loopout into a cell array of strings
    % ready to be written to a file
    
    cellLoops{1}=int2str(loopout.total.nloops);
    cellLoops{2}=int2str([loopout.total.nvertex, loopout.total.nfaces]);
    kk=3;
    
    % run through the different surfaces
    for ii=1:loopout.total.nloops
        cellLoops{kk}=int2str(loopout.surf(ii).nvertex);
        kk=kk+1;
        % Write an array of numbers in the correct format
        for jj=1:length(loopout.surf(ii).coord(:,1))
            cellLoops{kk}='';
            for ll=1:length(loopout.surf(ii).coord(jj,:))
                cellLoops{kk}=[cellLoops{kk},...
                    num2str(loopout.surf(ii).coord(jj,ll),24),'  '];
            end
            kk=kk+1;
        end
        
    end
    
end

function [isInternal]=FindInternalLoop(loop)
    % this function works by checkking all the inequalities around a
    % polygon
    isInternal=zeros(size(loop));
    for ii=1:length(loop)
        polygonPoints=loop(ii).subdivision;
%         polyVec=polygonPoints([2:end,1],:)-polygonPoints;
%         polyNorm=[polyVec(:,2),-polyVec(:,1)];
%         
%         numCond=size(polygonPoints,1);
        for jj=[1:ii-1,ii+1:length(loop)]
            
%             comparePoint=ones([size(polyNorm,1),1])*loop(jj).subdivision(1,:);
%             conditionNum=sum(polyNorm.*(comparePoint-polygonPoints),2);
%             allPos=sum(conditionNum>=0)==numCond;
%             allNeg=sum(conditionNum<=0)==numCond;
              isIn=inpolygon(loop(jj).subdivision(1,1),loop(jj).subdivision(1,2),...
                  loop(ii).subdivision(:,1),loop(ii).subdivision(:,2));
            isInternal(jj)=isInternal(jj) + (isIn);
        end
    end
    
end
%% Video Output functions

function []=MakeVideo(movStruct,fps,quality,fileName)
    
    writerObj = VideoWriter(fileName);
    writerObj.FrameRate=fps;
    writerObj.Quality=quality;
    open(writerObj)
    writeVideo(writerObj,movStruct)
    close(writerObj)
end

%% Layout Operation functions

function []=PersnaliseLayFile(FID,pltFile)
    
    frewind(FID);
    layData{1}='#!MC 1410';
    layData{2}=['$!VarSet |LFDSFN1| = ''"',pltFile,'"'''];
    if iscell(pltFile)
        for ii=1:length(pltFile)
            layData{1+ii}=['$!VarSet |LFDSFN',int2str(ii),'| = ''"',pltFile{ii},'"'''];
        end
    end
    WriteToFile(layData,FID)
end

%% File Opening Functions

function [FID]=OpenBoundaryFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['boundary_',marker,'.dat'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID,fileName]=OpenTecPLTFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['tec360dat_',marker,'.plt'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID]=OpenTecLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['tec360lay_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_Template.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfile(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function [FID]=OpenParamFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['param_',marker,'.dat'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID]=OpenCommentsFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['Comments_',marker,'.txt'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID]=OpenIndexFile(resultRoot,archiveName)
    % Creates a file in the current directory to write data to.
    
    resultRoot=MakePathCompliant(resultRoot);
    fileName=['Index_',archiveName,'.txt'];
    FID=fopen([resultRoot,filesep,archiveName,filesep,fileName],'a');
    
end

function [FID]=NameVideoFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['Video_',marker,'.avi'];
    FID=[writeDirectory,filesep,fileName];
    
end

%% Comments File

function [indexEntry]=MakeCommentsFile(FID,param,t,resultDirectory)
    
    varExtract={'typDat','case','noteFiles','tags'};
    [typDat,caseStr,noteFiles,tags]=ExtractVariables(varExtract,param);
    
    
    headerLines=GenerateCommentHeader(t,resultDirectory,typDat,caseStr,tags);
    automatedComments=ConcatenateAutomaticComments(noteFiles);
    indexEntry{1}=GenerateIndexEntry(t,resultDirectory,typDat,caseStr,tags);
    
    kk=1;
    breakLines{kk}=['--------------------------------------------------'];
    kk=kk+1;
    breakLines{kk}=['ADDITIONAL COMMENTS:'];
    kk=kk+1;
    breakLines{kk}=[' '];
    
    WriteToFile(headerLines,FID)
    WriteToFile(automatedComments,FID)
    WriteToFile(breakLines,FID)
    
end

function automatedComments=ConcatenateAutomaticComments(noteFiles)
    
    for ii=length(noteFiles):-1:1
        noteFileName{ii}=MakePathCompliant([cd,'\Result_Template\Notes_',noteFiles{ii},'.txt']);
        fidNote=fopen(noteFileName{ii},'r');
        rawComments(ii)=textscan(fidNote,'%s','Delimiter','\n');
        fclose(fidNote);
    end
    kk=1;
    automatedComments{kk}=['--------------------------------------------------'];
    kk=kk+1;
    automatedComments{kk}=['PRESET COMMENTS:'];
    kk=kk+1;
    automatedComments{kk}=[' '];
    for ii=1:length(noteFiles)
        kk=kk+1;
        automatedComments{kk}=['-------'];
        kk=kk+1;
        automatedComments{kk}=noteFileName{ii};
        kk=kk+1;
        automatedComments{kk}=['-------'];
        kk=kk+1;
        kk2=kk+length(rawComments{ii})-1;
        automatedComments(kk:kk2)=rawComments{ii};
        kk=kk2;
        kk=kk+1;
        automatedComments{kk}=[' '];
        
    end
    kk=kk+1;
    automatedComments{kk}=[' '];
    
end

function headerLines=GenerateCommentHeader(t,resultDirectory,typDat,caseStr,tags)
    
    kk=1;
    headerLines{kk}=datestr(t);
    kk=kk+1;
    headerLines{kk}=resultDirectory;
    kk=kk+1;
    headerLines{kk}=['Data File: ',typDat];
    kk=kk+1;
    headerLines{kk}=['Case: ',caseStr];
    kk=kk+1;
    headerLines{kk}=[' '];
    kk=kk+1;
    headerLines{kk}=['Tags: '];
    for ii=1:length(tags)
        headerLines{kk}=[headerLines{kk},tags{ii},', '];
    end
    kk=kk+1;
    headerLines{kk}=[' '];
    
    
end

function indexLine=GenerateIndexEntry(t,resultDirectory,typDat,caseStr,tags)
    
    indexLine='';
    indexLine=[indexLine,datestr(t)];
    indexLine=[indexLine,', '];
    indexLine=[indexLine,typDat];
    indexLine=[indexLine,', '];
    indexLine=[indexLine,caseStr];
    indexLine=[indexLine,', '];
    for ii=1:length(tags)-1
        indexLine=[indexLine,tags{ii},' - '];
    end
    indexLine=[indexLine,tags{ii+1}];
    indexLine=[indexLine,', '];
    indexLine=[indexLine,resultDirectory];
    
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

%% Error Report File

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

%% Generate Restart Binary

function []=GenerateRestartBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,filesep,'restart_',marker,'.mat'];
    save(fileName,'-struct','restartstruct');
    
end

function []=GenerateProfileBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,filesep,'restart_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'-struct','restartstruct');
    
end
