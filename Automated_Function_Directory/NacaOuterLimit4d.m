function [varargout]=NacaOuterLimit4d(varargin)
global NacaOuterLimit4d_Handle
nOut=nargout(NacaOuterLimit4d_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaOuterLimit4d_Handle(varargin{:});
end
