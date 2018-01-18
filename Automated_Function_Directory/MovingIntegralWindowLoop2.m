function [varargout]=MovingIntegralWindowLoop2(varargin)
% include_Utilities
global MovingIntegralWindowLoop2_Handle
nOut=nargout(MovingIntegralWindowLoop2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingIntegralWindowLoop2_Handle(varargin{:});
end
