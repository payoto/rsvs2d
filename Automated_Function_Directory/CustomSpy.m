function [varargout]=CustomSpy(varargin)
% include_PostProcessing
global CustomSpy_Handle
try
nOut=nargout(CustomSpy_Handle);
catch
include_PostProcessing
nOut=nargout(CustomSpy_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CustomSpy_Handle(varargin{:});
end
