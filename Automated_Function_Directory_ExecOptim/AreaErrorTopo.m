function [varargout]=AreaErrorTopo(varargin)
global AreaErrorTopo_Handle
nOut=nargout(AreaErrorTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=AreaErrorTopo_Handle(varargin{:});
end
