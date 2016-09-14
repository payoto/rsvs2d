%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a single snake validation case

function [procdat,out]=SnakeValid(validationName,parampreset,procdat)
    
    include_PostProcessing
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    include_PostProcessing
    include_Mex_Wrapper
    include_Validation
    
    runValid=true;
    if nargin==1
        parampreset=structInputVar('CurrentValidation');
        parampreset.structdat=GetStructureData(parampreset);
    end
    if nargin>2
        runValid=false;
    end
    
    param=ImposePresets(structInputVar('CurrentValidation'),parampreset);
    
    snakCell={'Snakestestsmooth1','Snakestestsmooth1_2','Snakestestsmooth2',...
        'Snakestestsmooth3','Snakestestsmooth3_1','Donught','Donught2',...
        'SnakesFoilVVSmall','BuzmanBiplane3','SnakesFoilVVSmall4',...
        'WeirdShapeIn','WeirdShapeOut'};
    
    T{length(snakCell)}=[];
    out=repmat(struct('unstructured',[],'loop',[],'unstructReshape',[],'snakSave',[]...
        ,'param',[],'rootdir',''),[1 length(snakCell)]);
    
    if runValid
        procdat=repmat(ProcDatTemplate,[1 length(snakCell)]);
        
        comStr=computer;
        
        if strcmp(comStr(1:2),'PC')
            [procdat,out,T]=WindowsSerialLoop(parampreset,snakCell,procdat,out);
        else
            [procdat,out,T]=LinuxParallelLoop(parampreset,snakCell,procdat,out);
        end
    end
        
    
    
    
    
    OutputDirectory(validationName,param,procdat,T)
    
    
end

function [procdat,out,T]=LinuxParallelLoop(parampreset,snakCell,procdat,out)
    
%     if numel(gcp('nocreate'))==0
%         poolName=parallel.importProfile('ExportOptimSnakesLinux.settings');
%         clusterObj=parcluster(poolName);
%         clusterObj.NumWorkers=12;
%         saveProfile(clusterObj);
%         parpool(poolName)
%     end
    for ii=1:length(snakCell)
        
        
        try 
            [T{ii},out(ii)]=CallMain(snakCell{ii},parampreset);
        catch ME
            T{ii}=ME.getReport;
        end
        [procdat(ii)]=ProcessData(out(ii),T{ii},snakCell{ii},procdat(ii));
    end
    
end

function [procdat,out,T]=WindowsSerialLoop(parampreset,snakCell,procdat,out)
    
    for ii=1:length(snakCell)
        try 
            [T{ii},out(ii)]=CallMain(snakCell{ii},parampreset);
            
        catch ME
            T{ii}=ME.getReport;
        end
        
        [procdat(ii)]=ProcessData(out(ii),T{ii},snakCell{ii},procdat(ii));
        
    end
    
end

function param=ImposePresets(param,parampreset)
    
    for ii=1:length(parampreset.structdat.vars)
        param=SetVariables({parampreset.structdat.vars(ii).name},...
            {ExtractVariables({parampreset.structdat.vars(ii).name},parampreset)},param);
    end
    
end

function [T,out]=CallMain(snakCell,parampreset)
    
    
    param=ImposePresets(structInputVar(['val_',snakCell]),parampreset);
    
    [T,out.unstructured,out.loop,out.unstructReshape,out.snakSave,out.param,out.rootdir]=...
        evalc('MainSimple(param)');
end

function [procdat]=ProcessData(out,T,caseName,procdat)
    % out is the output structure containing the data
    % T is the char array containing the text output
    
    lSnak=[out.snakSave(:).lSnak];
    DlSnak=lSnak(2:end)-lSnak(1:end-1);
    figure('Name',[caseName,'Length']), plot(1:length(lSnak),lSnak)
    figure('Name',[caseName,'LengthDelta']),
    semilogy(1:length(DlSnak),-(DlSnak),'o-',1:length(DlSnak),(DlSnak),'+-')
    
    try
        % read from text
        procdat.warningnum=numel(regexp(T,'Warning:'));
        
        procdat.errornum=numel(regexp(T,'Error '));
        
        timeMark=regexp(regexpi(T,'Iteration time[^\n]*','match'),':.*','match');
        procdat.cputime=str2num(regexprep(timeMark{1}{1},':',' '))*[3600;60;1;0.001]; % time between start and end. (look for "Iteration Time:")
        
        
        procdat.termination=~isempty(regexpi(T,'Snakes Converged!')); % true / false (true means terminate with convergence) (look for Snakes Converged!)
        
        % read from array
        procdat.caseName=caseName;
        procdat.niter=length(out.snakSave); % number of iteration to termination
        
        procdat.snaketime=sum([out.snakSave(:).dt]);
        procdat.volerror=log10(out.snakSave(end).currentConvVolume);
        procdat.velerror=log10(out.snakSave(end).currentConvVelocity);
        procdat.length=out.snakSave(end).lSnak;
        procdat.path=out.rootdir;
    catch
        procdat.errornum=1;
    end
    
end

function [tableCell]=GenerateTableResult(validationName,t,marker,procdat)
    
    paramCell{1}=['# ',validationName,' - ',datestr(t),' - ',marker];
    paramCell{2}=['# '];
    
    
    fieldsProc=fieldnames(procdat);
    nProc=length(fieldsProc);
    entryCell{1+length(procdat)}={''};
    
    entryCell{1}=fieldsProc{1};
    for ii=2:nProc
        entryCell{1}=[entryCell{1},' , ',fieldsProc{ii}];
    end
    
    for ii=1:length(procdat)
        str=ProcesstoString(procdat(ii).(fieldsProc{1}));
        for jj=2:nProc
            str=[str,' , ',ProcesstoString(procdat(ii).(fieldsProc{jj}))];
        end
        
        
        entryCell{ii+1}=str;
    end
    tableCell=[paramCell,entryCell];
    
end

function [sumLines]=GenerateSummaryEntry(validationName,t,marker,procdat)
    
    sumLines{1}=['# ',validationName,' - ',datestr(t),' - ',marker];
    
    [sumLines(2:3)]=GenerateSummaryLine(procdat);
    
    
end

function []=OutputDirectory(validationName,param,procdat,T)
    
    varExtract={'resultRoot','archiveName'};
    [resultRoot,archiveName]=ExtractVariables(varExtract,param);
    
    [marker,t]=GenerateResultMarker(validationName);
    [writeDirectory]=GenerateValidationDirectory(marker,resultRoot,archiveName,t);
    
    fileName=[writeDirectory,filesep,'Procdat_',marker,'.mat'];
    save(fileName,'procdat');
    
    fileName=[writeDirectory,filesep,'Diary_',marker,'.txt'];
    FIDDiary=fopen(fileName,'w');
    WriteToFile(T,FIDDiary);
    fclose(FIDDiary);
    
    % Parameter Data
    [fidParam]=OpenParamFile(writeDirectory,marker);
    GenerateParameterFile(fidParam,param,t,marker);
    
    [tableCell]=GenerateTableResult(validationName,t,marker,procdat);
    OutputTableResult(tableCell,writeDirectory,marker);
    
    [sumLines]=GenerateSummaryEntry(validationName,t,marker,procdat);
    OutputValidationSummary(resultRoot,archiveName,sumLines)
    
    OutputAllFigures(writeDirectory)
end

function []=OutputTableResult(tableCell,writeDirectory,marker)
    
    fileName=['table_',marker,'.csv'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
    WriteToFile(tableCell,FID)
    
    fclose(FID);
end

function []=OutputAllFigures(writeDirectory)
    h=findobj('type','figure');
    
    figDir=[writeDirectory,filesep,'fig'];
    mkdir(figDir)
    
    for ii=1:length(h)
        figName=[figDir,filesep,'fig',int2str(ii),'_',h(ii).Name,'.fig'];
        
        hgsave(h(ii),figName);
        
        
    end
    
    close all
end

function []=OutputValidationSummary(resultRoot,archiveName,sumLines)
    % Creates a file in the current directory to write data to.
    
    resultRoot=MakePathCompliant(resultRoot);
    fileName=['Validations_',archiveName,'.csv'];
    FID=fopen([resultRoot,filesep,archiveName,filesep,fileName],'a');
    
    WriteToFile(sumLines,FID)
    
    fclose(FID);
    
end

function [resultDirectory]=GenerateValidationDirectory(marker,resultRoot,...
        archiveName,t)
    if ~exist('t','var'),t=now;end
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,filesep,archiveName,filesep,dateSubFolders,...
        filesep,'VALIDATION_',marker];
    
    resultDirectory=MakePathCompliant(resultDirectory);
    
    mkdir(resultDirectory)
end

function [outLine]=GenerateValidationData(procdat)
    
    kk=1;
    outLine{kk,1}=sum([procdat(:).warningnum]);outLine{kk,2}='int';outLine{kk,3}='warnings';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).errornum]);outLine{kk,2}='int';outLine{kk,3}='errors';
    kk=kk+1;
    outLine{kk,1}=[sum([procdat(:).termination]),length(procdat)];
    outLine{kk,2}='int';outLine{kk,3}='termination';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).cputime]);outLine{kk,2}='double';outLine{kk,3}='total cputime';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).snaketime]);outLine{kk,2}='double';outLine{kk,3}='total snaketime';
    kk=kk+1;
    outLine{kk,1}=mean([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - Mean';
    kk=kk+1;
    outLine{kk,1}=std([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - STD';
    kk=kk+1;
    outLine{kk,1}=max([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - max';
    kk=kk+1;
    outLine{kk,1}=min([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - min';
    kk=kk+1;
    outLine{kk,1}=mean([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - Mean';
    kk=kk+1;
    outLine{kk,1}=std([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - STD';
    kk=kk+1;
    outLine{kk,1}=max([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - max';
    kk=kk+1;
    outLine{kk,1}=min([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - min';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).length]);outLine{kk,2}='double';outLine{kk,3}='total length';
    
end


%% Simplified Main File

function [unstructured,loop,unstructReshape,snakSave,param,rootDirectory]=MainSimple(param)
    % Main function for the execution of the Subdivision process
    
    
    diaryFile=[cd,'\Result_Template\Latest_Diary.log'];
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    
    [param,unstructured,unstructuredrefined,loop,connectstructinfo...
        ,snakSave,unstructReshape,gridrefined,restartsnake]=StandardRun(param);
    
    varExtract={'typeBound','refineSteps','case'};
    [typeBound,refineSteps,caseString]=ExtractVariables(varExtract,param);
    
    % Restart structure generation
    restartstruct.gridrefined=gridrefined;
    restartstruct.param=param;
    restartstruct.connectstructinfo=connectstructinfo;
    restartstruct.unstructReshape=unstructReshape;
    restartstruct.snakrestart=restartsnake;
    
    % Post processes
    loop=SubdivisionSurface_Snakes(loop,refineSteps,param);
    CheckResults(unstructured,loop,typeBound,caseString)
    
    tecoutstruct.baseGrid=unstructReshape;
    tecoutstruct.fineGrid=unstructuredrefined;
    tecoutstruct.snakSave=snakSave;
    tecoutstruct.connectstructinfo=connectstructinfo;
    
    diary off
    rootDirectory=ManageOutputResults(param,loop,tecoutstruct,restartstruct);
    %TecplotOutput(unstructReshape,unstructuredrefined,snakSave,connectstructinfo)
    %OutPutBinaryResults(snakSave,saveParam,typDat)
    
    
end

%% Top Level Execution processes

function [param,unstructured,unstructuredrefined,loop,connectstructinfo...
        ,snakSave,unstructReshape,gridrefined,restartsnake]...
        =StandardRun(param)
    
    
    param.general.restart=false;
    varExtract={'useSnakes'};
    [useSnakes]=ExtractVariables(varExtract,param);
    
    
    % Execution of subroutines
    
    [unstructured,loop,unstructReshape]=ExecuteGridInitialisation(param);
    
    [gridrefined,connectstructinfo,unstructuredrefined,looprefined]=...
        ExecuteGridRefinement(unstructReshape,param);
    snakSave=[];
    if useSnakes
        [snaxel,snakposition,snakSave,loop,restartsnake]=ExecuteSnakes(gridrefined,looprefined,...
            unstructReshape,connectstructinfo,param);
    end
    
end

function [unstructured,loop,unstructReshape]...
        =ExecuteGridInitialisation(param)
    % Executes the Grid Initialisation process
    
    disp('GENERATION PROCESS START')
    t1=now;
    [unstructured,loop,unstructReshape]=...
        GridInitialisationV2(param);
    disp('Grid Initialisation exited')
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GENERATION PROCESS end')
    
    
end

function [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        ExecuteGridRefinement(unstructReshape,param)
    % Executes the Grid Initialisation process
    
    disp('GRID REFINEMENT START')
    t1=now;
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GRID REFINEMENT end')
    
    
end

function [snaxel,snakposition,snakSave,loop,restartsnake]=ExecuteSnakes(unstructured,loop,...
        oldGrid,connectionInfo,param)
    % Executes the snakes edge detection process
    %
    
    t1=now;
    disp('SNAKE PROCESS START')
    [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=Snakes(unstructured,loop,...
        oldGrid,connectionInfo,param);
    
    t2=now;
    disp(['Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('SNAKE PROCESS END')
    
    if length(loopsnaxel)==length(loop)
        for ii=1:length(loopsnaxel)
            loop(ii).snaxel=loopsnaxel(ii).snaxel;
        end
    else
        loop=loopsnaxel;
    end
end


%% Subdivision process

function [loop]=SubdivisionSurface_Snakes(loop,refineSteps,param)
    % Function taking in a closed loop of vertices and applying the subdivision
    % process
    % typBOund is te type of boundary that is expected, it can either be the
    % string 'vertex' (default) or the string 'snaxel' to show that the snaxel
    % process has been used
    varExtract={'typeBound','subdivType','TEShrink','LEShrink','typeCorner'};
    [typeBound,subdivType,TEShrink,LEShrink,typeCorner]=ExtractVariables(varExtract,param);
    if ~exist('typeBound','var'), typeBound='vertex'; end
    
    sharpen=[LEShrink,TEShrink];
    
    for ii=1:length(loop)
        startPoints=loop(ii).(typeBound).coord;
        loop(ii).isccw=CCWLoop(startPoints);
        newPoints=SubDivision(startPoints,refineSteps,subdivType,sharpen,typeCorner);
        %newPoints=SubSurfBSpline(startPoints(1:end-2,:),refineSteps);
        loop(ii).subdivision=newPoints;
    end
end


%% Plot Functions
function []=CheckResults(unstructured,loop,typeBound,caseString)
    global nDim domainBounds
    
    if nDim==2
        figh=figure('Name',['Final Profile ',caseString]);
        axh=axes;
        hold on
        
        colString='bgcmyk';
        
        isEdgeIndex=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeIndex)
            PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'bo')
        end
        
        isEdgeIndex=find(unstructured.edge.boundaryis0);
        for ii=1:length(isEdgeIndex)
            PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'b-')
        end
        
        
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        
        for ii=1:length(loop)
            [~,colIndex]=IntegerQuotient(ii,length(colString));
            colIndex=colIndex+1;
            PlotLoop(figh,axh,loop,ii,[colString(colIndex),'--'],typeBound)
            PlotSubDiv(figh,axh,loop,ii,[colString(colIndex),'-'])
        end
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    end
    
end

function []=PlotEdge(figh,axh,unstructured,indexEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(indexEdge,:);
    coord=unstructured.vertex.coord(vertices,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotLoop(figh,axh,loop,indexLoop,format,typeBound)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).(typeBound).coord;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotSubDiv(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).subdivision;
    
    plot(coord(:,1),coord(:,2),format)
    
end
