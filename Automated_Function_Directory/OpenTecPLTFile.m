function [varargout]=OpenTecPLTFile(varargin)
% include_PostProcessing
global OpenTecPLTFile_Handle
try
nOut=nargout(OpenTecPLTFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenTecPLTFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenTecPLTFile_Handle(varargin{:});
end
