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
    try 
        marker=param.outstart.marker;
        t=param.outstart.t;
        rootDirectory=param.outstart.rootDirectory;
        writeDirectory=rootDirectory;
    catch
        
    [marker,t]=GenerateResultMarker(matlab.lang.makeValidName(caseStr));
    [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t);
    rootDirectory=writeDirectory;
    end
    
    % Create Marker
    % Create Directory
    CopyDiary(writeDirectory,marker)
    GenerateRestartBinary(writeDirectory,marker,restartstruct)
    
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
    
    TecplotOutput('snakes',fidTecPLT,baseGrid,fineGrid,snakSave,connectstructinfo)
    
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
    
    
    fclose('all');
end

function [rootDirectory,param]=ManageResultsStart(param)
    
    varExtract={'makeMov','typDat','resultRoot','archiveName','case'};
    [makeMov,typDat,resultRoot,archiveName,caseStr]=ExtractVariables(varExtract,param);
    [param.outstart.marker,param.outstart.t]=GenerateResultMarker(matlab.lang.makeValidName(caseStr));
    [rootDirectory]=GenerateResultDirectoryName(param.outstart.marker,resultRoot,archiveName,param.outstart.t);
    param.outstart.rootDirectory=rootDirectory;
end


