%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%          Function for Creation of Header
%                     Folders
%                Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [] = HeaderActivation(funcHandles,funcDir)
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    c=evalc('display(funcHandles)');
    pattern='@\w*';
    
    funcHandlesNamesCell=regexp(c,pattern,'match');
    system(['mkdir "',funcDir,'"']);
    for ii=1:length(funcHandlesNamesCell)
        funcName=funcHandlesNamesCell{ii}(2:end);
        funcHandleVarname=[funcHandlesNamesCell{ii}(2:end),'_Handle'];
        eval(['global ',funcHandleVarname])
        eval([funcHandleVarname,'=funcHandles{ii};']);
        WriteConatainerFunctionFile(funcName,funcDir);
    end
    addpath(funcDir);
end

function []=WriteConatainerFunctionFile(funcName,funcDir)
    
    fID=fopen([funcDir,'\',funcName,'.m'],'w');
    
    funcText{1}=['function [varargout]=',funcName,'(varargin)'];
    funcText{2}=['global ',funcName,'_Handle'];
    funcText{3}=['nOut=nargout(',funcName,'_Handle',');'];
    funcText{4}=['nOutReq=nargout;'];
    funcText{5}=['nOut(nOut<0)=nOutReq;'];
    funcText{6}=['[varargout{1:nOut}]=',funcName,'_Handle','(varargin{:});'];
    funcText{7}=['end'];
    
    for ii=1:length(funcText)
        
        fprintf(fID,[funcText{ii},'\n']);
        
    end
    fclose('all');
    
end


