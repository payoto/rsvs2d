function [varargout]=OpenParamFile(varargin)
% include_PostProcessing
global OpenParamFile_Handle
try
nOut=nargout(OpenParamFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenParamFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenParamFile_Handle(varargin{:});
end
