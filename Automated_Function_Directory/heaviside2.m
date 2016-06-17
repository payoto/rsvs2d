function [varargout]=heaviside2(varargin)
global heaviside2_Handle
nOut=nargout(heaviside2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=heaviside2_Handle(varargin{:});
end
