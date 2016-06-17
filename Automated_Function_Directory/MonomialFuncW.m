function [varargout]=MonomialFuncW(varargin)
global MonomialFuncW_Handle
nOut=nargout(MonomialFuncW_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MonomialFuncW_Handle(varargin{:});
end
