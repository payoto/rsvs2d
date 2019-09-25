function [varargout]=OpenCommentsFile(varargin)
% include_PostProcessing
global OpenCommentsFile_Handle
try
nOut=nargout(OpenCommentsFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenCommentsFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenCommentsFile_Handle(varargin{:});
end
