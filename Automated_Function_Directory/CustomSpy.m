function [varargout]=CustomSpy(varargin)
% include_PostProcessing
global CustomSpy_Handle
nOut=nargout(CustomSpy_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CustomSpy_Handle(varargin{:});
end
