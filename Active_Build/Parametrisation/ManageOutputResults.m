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



function [rootDirectory,varargout]=ManageOutputResults(entryPoint,param,loop,tecoutstruct,restartstruct)
    
    if ~exist('entryPoint','var');
        entryPoint='end';
    end
    varargout{1}=param;
    
    switch entryPoint
        
        case 'end'
            [rootDirectory]=ManageResultsEnd(param,loop,tecoutstruct,restartstruct);
        case 'start'
            [rootDirectory,param]=ManageResultsStart(param);
    end
    varargout{1}=param;
    
    
    
end

function [rootDirectory]=ManageResultsEnd(param,loop,tecoutstruct,restartstruct)
    % Unpack necessary variables
    varExtract={'makeMov','typDat','resultRoot','archiveName','case'};
    [makeMov,typDat,resultRoot,archiveName,caseStr]=ExtractVariables(varExtract,param);
    kk=0;
    try
        marker=param.outstart.marker;
        t=param.outstart.t;
        rootDirectory=param.outstart.rootDirectory;
        writeDirectory=rootDirectory;
    catch
        
        [marker,t]=GenerateResultMarker(matlab.lang.makeValidName(caseStr));
        [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t);
        rootDirectory=writeDirectory;
        try
            CopyDiary(writeDirectory,marker)
        catch ME
            kk=kk+1;
            T{kk}=ME;
        end
    end
    
    % Create Marker
    % Create Directory
    
    
    try
        GenerateRestartBinary(writeDirectory,marker,restartstruct)
    catch ME
        kk=kk+1;
        T{kk}=ME;
    end
    % Output boundary data file
    try
        [fidBoundary]=OpenBoundaryFile(writeDirectory,marker);
        BoundaryOutput(loop,fidBoundary);
        fclose(fidBoundary);
    catch ME
        kk=kk+1;
        T{kk}=ME;
    end
    % TecPlot Data
    try
        [fidTecPLT,pltFileName]=OpenTecPLTFile(writeDirectory,marker);
        fidTecLAY=OpenTecLayFile(writeDirectory,marker);
        PersnaliseLayFile(fidTecLAY,pltFileName)
        
        baseGrid=tecoutstruct.baseGrid;
        fineGrid=tecoutstruct.fineGrid;
        snakSave=tecoutstruct.snakSave;
        connectstructinfo=tecoutstruct.connectstructinfo;
        
        TecplotOutput('snakes',fidTecPLT,baseGrid,fineGrid,snakSave,connectstructinfo)
    catch ME
        kk=kk+1;
        T{kk}=ME;
    end
    % Parameter Data
    [fidParam]=OpenParamFile(writeDirectory,marker);
    GenerateParameterFile(fidParam,param,t,marker);
    
    
    % Comments file
    try
        [fidComments]=OpenCommentsFile(writeDirectory,marker);
        indexEntry=MakeCommentsFile(fidComments,param,t,writeDirectory);
        fclose(fidComments);
    catch ME
        kk=kk+1;
        T{kk}=ME;
    end
    % Video File
    try
        if makeMov
            fps=5;
            quality=100;
            movStruct=[snakSave(:).movFrame];
            [fileName]=NameVideoFile(writeDirectory,marker);
            MakeVideo(movStruct,fps,quality,fileName);
        end
    catch ME
        kk=kk+1;
        T{kk}=ME;
    end
    % Index File Entry
    try
        [fidIndexFile]=OpenIndexFile(resultRoot,archiveName);
        WriteToFile(indexEntry,fidIndexFile);
        fclose(fidIndexFile);
    catch ME
        kk=kk+1;
        T{kk}=ME;
    end
    % Generate Restart Binary
    
    
    fclose('all');
    
    if kk>0
        for ii=1:kk
            disp('')
            disp('------------------------------')
            disp('------------------------------')
            disp(T{ii}.getReport)
        end
        warning('Snake Data Output was partially succesful')
    end
    
    
end

function [rootDirectory,outstart]=ManageResultsStart(param)
    
    varExtract={'makeMov','typDat','resultRoot','archiveName','case'};
    [makeMov,typDat,resultRoot,archiveName,caseStr]=ExtractVariables(varExtract,param);
    [outstart.marker,outstart.t]=GenerateResultMarker(matlab.lang.makeValidName(caseStr));
    [rootDirectory]=GenerateResultDirectoryName(outstart.marker,resultRoot,archiveName,outstart.t);
    outstart.rootDirectory=rootDirectory;
end


