function [varargout]=CompareProfilesAreaSquaredTopo(varargin)
% InverseDesign_ErrorTopo
global CompareProfilesAreaSquaredTopo_Handle
nOut=nargout(CompareProfilesAreaSquaredTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesAreaSquaredTopo_Handle(varargin{:});
end
