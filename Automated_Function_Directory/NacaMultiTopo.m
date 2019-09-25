function [varargout]=NacaMultiTopo(varargin)
% include_Optimisation
global NacaMultiTopo_Handle
try
nOut=nargout(NacaMultiTopo_Handle);
catch
include_Optimisation
nOut=nargout(NacaMultiTopo_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaMultiTopo_Handle(varargin{:});
end
