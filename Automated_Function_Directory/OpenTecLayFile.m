function [varargout]=OpenTecLayFile(varargin)
% include_PostProcessing
global OpenTecLayFile_Handle
try
nOut=nargout(OpenTecLayFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenTecLayFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenTecLayFile_Handle(varargin{:});
end
