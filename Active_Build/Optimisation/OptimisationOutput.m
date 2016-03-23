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


function []=OptimisationOutput(paramoptim)
    
    % Unpack necessary variables
    varExtract={'optimCase','typDat','resultRoot','archiveName'};
    [optimCase,typDat,resultRoot,archiveName]=ExtractVariables(varExtract,paramoptim);
    
    [marker,t]=GenerateResultMarker([optimCase,'_',typDat]);
    [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t);
    
    CopyDiary(writeDirectory,marker)
    
end


%% 