function [procdat,out]=SnakeValid(validationName,procdat,out)
    
    include_PostProcessing
    
    snakCell={'Snakestestsmooth1','Snakestestsmooth1_2','Snakestestsmooth2',...
        'Snakestestsmooth3','Snakestestsmooth3_1','Donught','Donught2',...
        'SnakesFoilVVSmall','BuzmanBiplane3','SnakesFoilVVSmall4',...
        'WeirdShapeIn','WeirdShapeOut'};
    
    T{length(snakCell)}=[];
    if nargin==1
        out=repmat(struct('unstructured',[],'loop',[],'unstructReshape',[],'snakSave',[]...
            ,'param',[],'rootdir',''),[1 length(snakCell)]);
        procdat=repmat(struct('caseName','','warningnum',0,'errornum',0,'cputime',0,'termination',false, 'niter',0,...
            'snaketime',0,'volerror',1,'velerror',1,'length',0,'path',''),[1 length(snakCell)]);
        
        for ii=1:length(snakCell)
            pause(ii/5)
            [T{ii},out(ii)]=CallMain(snakCell{ii});
            [procdat(ii)]=ProcessData(out(ii),T{ii},snakCell{ii},procdat(ii));
        end
    end
    param=structInputVar('CurrentValidation');
    
    OutputDirectory(validationName,param,procdat)
    
    
end

function [T,out]=CallMain(snakCell)
    
    [T,out.unstructured,out.loop,out.unstructReshape,out.snakSave,out.param,out.rootdir]=...
        evalc('Main([''val_'',snakCell])');
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
        procdat.warningnum=numel(regexpi(T,'warning'));
        
        procdat.errornum=numel(regexpi(T,'error'));
        
        timeMark=regexp(regexpi(T,'Iteration time[^\n]*','match'),':.*','match');
        procdat.cputime=str2num(regexprep(timeMark{1}{1},':',' '))*[3600;60;1;0.001]; % time between start and end. (look for "Iteration Time:")
        
        procdat.termination=~isempty(regexpi(T,'Snakes Converged!')); % true / false (true means terminate with convergence) (look for Snakes Converged!)
        
        % read from array
        procdat.caseName=caseName;
        procdat.niter=length(out.snakSave); % number of iteration to termination
        
        procdat.snaketime=sum([out.snakSave(:).dt]);
        procdat.volerror=out.snakSave(end).currentConvVolume;
        procdat.velerror=out.snakSave(end).currentConvVelocity;
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

function [strOut]=ProcesstoString(inputVar,classin)
    
    if nargin==1
        classin=class(inputVar);
    end
    
    switch classin
        case 'int'
            strOut=int2str(inputVar);
        case 'logical'
            strOut=int2str(inputVar);
        case 'double'
            if all(mod(inputVar,1)==0)
                strOut=int2str(inputVar);
            else
                strOut=num2str(inputVar,' %12.7e ');
            end
        case 'char'
            strOut=inputVar;
        otherwise
            strOut='';
            
    end
    
    
    
    
end

function [sumLines]=GenerateSummaryLine(validationName,t,marker,procdat)
    
    [outLine]=GenerateOutputData(procdat);
    
    
    sumLines{1}=['# ',validationName,' - ',datestr(t),' - ',marker];
    
    separator='';
    sumLines{2}='# ';
    sumLines{3}='';
    for ii=1:size(outLine,1)
        sumLines{2}=[sumLines{2},separator,outLine{ii,3}];
        sumLines{3}=[sumLines{3},separator,ProcesstoString(outLine{ii,1},outLine{ii,2})];
        separator=' , ';
    end
    

end

function []=OutputDirectory(validationName,param,procdat)
    
    varExtract={'resultRoot','archiveName'};
    [resultRoot,archiveName]=ExtractVariables(varExtract,param);
    
    
    [marker,t]=GenerateResultMarker(validationName);
    [writeDirectory]=GenerateValidationDirectory(marker,resultRoot,archiveName,t);
    
    % Parameter Data
    [fidParam]=OpenParamFile(writeDirectory,marker);
    GenerateParameterFile(fidParam,param,t,marker);
    
    [tableCell]=GenerateTableResult(validationName,t,marker,procdat);
    OutputTableResult(tableCell,writeDirectory,marker);
    
    [sumLines]=GenerateSummaryLine(validationName,t,marker,procdat);
    OutputValidationSummary(resultRoot,archiveName,sumLines)
end

function []=OutputTableResult(tableCell,writeDirectory,marker)
    
    fileName=['table_',marker,'.dat'];
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
    
end

function []=OutputValidationSummary(resultRoot,archiveName,sumLines)
    % Creates a file in the current directory to write data to.
    
    resultRoot=MakePathCompliant(resultRoot);
    fileName=['Validations_',archiveName,'.txt'];
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

function [outLine]=GenerateOutputData(procdat)
    
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