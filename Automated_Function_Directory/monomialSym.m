function [varargout]=monomialSym(varargin)
global monomialSym_Handle
nOut=nargout(monomialSym_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=monomialSym_Handle(varargin{:});
end
