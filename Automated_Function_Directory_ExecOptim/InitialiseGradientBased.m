function [varargout]=InitialiseGradientBased(varargin)
global InitialiseGradientBased_Handle
nOut=nargout(InitialiseGradientBased_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseGradientBased_Handle(varargin{:});
end
