function [varargout]=NacaOuterLimit(varargin)
global NacaOuterLimit_Handle
nOut=nargout(NacaOuterLimit_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaOuterLimit_Handle(varargin{:});
end
