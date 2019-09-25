function [varargout]=OpenIndexFile(varargin)
% include_PostProcessing
global OpenIndexFile_Handle
try
nOut=nargout(OpenIndexFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenIndexFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenIndexFile_Handle(varargin{:});
end
