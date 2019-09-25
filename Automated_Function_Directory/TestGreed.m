function [varargout]=TestGreed(varargin)
% include_PostProcessing
global TestGreed_Handle
try
nOut=nargout(TestGreed_Handle);
catch
include_PostProcessing
nOut=nargout(TestGreed_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestGreed_Handle(varargin{:});
end
