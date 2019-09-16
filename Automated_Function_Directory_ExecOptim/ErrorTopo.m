function [varargout]=ErrorTopo(varargin)
% InverseDesign_ErrorTopo
global ErrorTopo_Handle
nOut=nargout(ErrorTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ErrorTopo_Handle(varargin{:});
end
