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
    disp('WRITING AUTOMATIC FUNCTIONS')
    funcHandlesNamesCell=regexp(c,pattern,'match');
    funcDir=MakePathCompliant(funcDir);
    if ~isdir(funcDir)
        
        mkdir(funcDir);
    end
    for ii=1:length(funcHandlesNamesCell)
        funcName=funcHandlesNamesCell{ii}(2:end);
        funcHandleVarname=[funcHandlesNamesCell{ii}(2:end),'_Handle'];
        eval(['global ',funcHandleVarname])
        eval([funcHandleVarname,'=funcHandles{ii};']);
        WriteContainerFunctionFile(funcName,funcDir);
    end
    addpath(funcDir);
end

function pathName=MakePathCompliant(pathName)
    
    compStr=computer;
    if strcmp(compStr(1:2),'PC')
        pathName=regexprep(pathName,'/','\\');
    else
        
        pathName=regexprep(pathName,'\\','/');
    end
end

function []=WriteContainerFunctionFile(funcName,funcDir)
    fID=-1;
    t=now;
    Dt=0;
    while(fID<0 && Dt<0.01)
        fID=fopen(MakePathCompliant([funcDir,'\',funcName,'.m']),'w');
        Dt=now-t;
    end
    
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


