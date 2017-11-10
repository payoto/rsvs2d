function [varargout]=LSMPoly(varargin)
global LSMPoly_Handle
nOut=nargout(LSMPoly_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LSMPoly_Handle(varargin{:});
end
