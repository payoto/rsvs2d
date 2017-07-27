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



function [rootDirectory]=ManageOutputResults(param,loop,tecoutstruct,restartstruct)
    
    % Unpack necessary variables
    varExtract={'makeMov','typDat','resultRoot','archiveName','case'};
    [makeMov,typDat,resultRoot,archiveName,caseStr]=ExtractVariables(varExtract,param);
    
    
    % Create Marker
    [marker,t]=GenerateResultMarker(matlab.lang.makeValidName(caseStr));
    % Create Directory
    [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t);
    rootDirectory=writeDirectory;
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
    GenerateRestartBinary(writeDirectory,marker,restartstruct)
    
    
    fclose('all');
    
end
