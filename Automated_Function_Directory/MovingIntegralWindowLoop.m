function [varargout]=MovingIntegralWindowLoop(varargin)
% include_Utilities
global MovingIntegralWindowLoop_Handle
nOut=nargout(MovingIntegralWindowLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingIntegralWindowLoop_Handle(varargin{:});
end
