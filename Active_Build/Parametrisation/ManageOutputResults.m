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



function []=ManageOutputResults(param,loop,tecoutstruct,restartstruct)
    
    % Unpack necessary variables
    varExtract={'makeMov','typDat','resultRoot','archiveName'};
    [makeMov,typDat,resultRoot,archiveName]=ExtractVariables(varExtract,param);
    
    
    % Create Marker
    [marker,t]=GenerateResultMarker(typDat);
    % Create Directory
    [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t);
    
    CopyDiary(writeDirectory,marker)
    
    % Output boundary data file
    [fidBoundary]=OpenBoundaryFile(writeDirectory,marker);
    BoundaryOutput(loop,fidBoundary);
    fclose(fidBoundary);
    
    % TecPlot Data
    [fidTecPLT,pltFileName]=OpenTecPLTFile(writeDirectory,marker);
    fidTecLAY=OpenTecLayFile(writeDirectory,marker);
    PersnaliseLayFile(fidTecLAY,pltFileName)
    
    baseGrid=tecoutstruct.baseGrid;
    fineGrid=tecoutstruct.fineGrid;
    snakSave=tecoutstruct.snakSave;
    connectstructinfo=tecoutstruct.connectstructinfo;
    
    TecplotOutput(fidTecPLT,baseGrid,fineGrid,snakSave,connectstructinfo)
    
    % Parameter Data
    [fidParam]=OpenParamFile(writeDirectory,marker);
    GenerateParameterFile(fidParam,param,t,marker);
    
    
    % Comments file
    [fidComments]=OpenCommentsFile(writeDirectory,marker);
    indexEntry=MakeCommentsFile(fidComments,param,t,writeDirectory);
    fclose(fidComments);
    % Video File
    
    if makeMov
        fps=5;
        quality=100;
        movStruct=[snakSave(:).movFrame];
        [fileName]=NameVideoFile(writeDirectory,marker);
        MakeVideo(movStruct,fps,quality,fileName);
    end
    
    % Index File Entry 
    [fidIndexFile]=OpenIndexFile(resultRoot,archiveName);
    WriteToFile(indexEntry,fidIndexFile);
    fclose(fidIndexFile);
    % Generate Restart Binary
    GenerateRestartBinary(writeDirectory,marker,restartstruct)
    
    
    fclose('all');
    
end


%% Boundary Output to .dat file

function []=BoundaryOutput(loop,FID)
    
    % trim loops and extract data
    loopout=TrimLoops(loop);
    % format numeric data to printable strings
    cellLoops=DataToString(loopout);
    
    % print string data to file
    WriteToFile(cellLoops,FID)
    
end

function [loopout]=TrimLoops(loop)
    % function extracting the data that must be written to the boundary.dat
    % file
    
    nLoop=length(loop);
    
    for ii=1:nLoop
        if loop(ii).isccw
            loopout.surf(ii).coord=loop(ii).subdivision(1:end-2,:);
        else
            loopout.surf(ii).coord=loop(ii).subdivision(end-2:-1:1,:);
        end
        loopout.surf(ii).nvertex=length(loopout.surf(ii).coord(:,1));
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
                    num2str(loopout.surf(ii).coord(jj,ll),8),'  '];
            end
            kk=kk+1;
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
    WriteToFile(layData,FID)
end


%% File Opening Functions

function [FID]=OpenBoundaryFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['boundary_',marker,'.dat'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID,fileName]=OpenTecPLTFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['tec360dat_',marker,'.plt'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenTecLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['tec360lay_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_Template.lay'];
    copyfile(originalLayFile,[writeDirectory,'\',fileName])
    FID=fopen([writeDirectory,'\',fileName],'r+');
    
end

function []=CopyDiary(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['diary_',marker,'.log'];
    originalLayFile=[cd,'\Result_Template\Latest_Diary.log'];
    copyfile(originalLayFile,[writeDirectory,'\',fileName])
    
end

function [FID]=OpenParamFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['param_',marker,'.dat'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenCommentsFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['Comments_',marker,'.txt'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenIndexFile(resultRoot,archiveName)
    % Creates a file in the current directory to write data to.
    
    fileName=['Index_',archiveName,'.txt'];
    FID=fopen([resultRoot,'\',archiveName,'\',fileName],'a');
    
end

function [FID]=NameVideoFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['Video_',marker,'.avi'];
    FID=[writeDirectory,'\',fileName];
    
end

%% Open Result Directory

% Now global functions

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
        noteFileName{ii}=[cd,'\Result_Template\Notes_',noteFiles{ii},'.txt'];
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

%% parameter File

% Now Global Functions

%% Generate Restart Binary

function []=GenerateRestartBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,'\restart_',marker,'.mat'];
    save(fileName,'-struct','restartstruct');
    
end

