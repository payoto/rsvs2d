function [varargout]=v(varargin)
global v_Handle
nOut=nargout(v_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=v_Handle(varargin{:});
end
