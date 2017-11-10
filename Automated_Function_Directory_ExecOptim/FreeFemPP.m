function [varargout]=FreeFemPP(varargin)
global FreeFemPP_Handle
nOut=nargout(FreeFemPP_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FreeFemPP_Handle(varargin{:});
end
