function [varargout]=MovingIntegralWindowLoop(varargin)
% include_Utilities
global MovingIntegralWindowLoop_Handle
try
nOut=nargout(MovingIntegralWindowLoop_Handle);
catch
include_Utilities
nOut=nargout(MovingIntegralWindowLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingIntegralWindowLoop_Handle(varargin{:});
end
