function [varargout]=TestGreed(varargin)
% include_PostProcessing
global TestGreed_Handle
nOut=nargout(TestGreed_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestGreed_Handle(varargin{:});
end
