function [varargout]=NacaMultiTopo(varargin)
global NacaMultiTopo_Handle
nOut=nargout(NacaMultiTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaMultiTopo_Handle(varargin{:});
end
