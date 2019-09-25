function [varargout]=MovingIntegralWindowLoop2(varargin)
% include_Utilities
global MovingIntegralWindowLoop2_Handle
try
nOut=nargout(MovingIntegralWindowLoop2_Handle);
catch
include_Utilities
nOut=nargout(MovingIntegralWindowLoop2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingIntegralWindowLoop2_Handle(varargin{:});
end
