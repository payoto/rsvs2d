function [varargout]=GenerateResultDirectoryName(varargin)
% include_PostProcessing
global GenerateResultDirectoryName_Handle
try
nOut=nargout(GenerateResultDirectoryName_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateResultDirectoryName_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateResultDirectoryName_Handle(varargin{:});
end
