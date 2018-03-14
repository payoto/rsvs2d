function [varargout]=CompareProfilesAreaTopo(varargin)
% InverseDesign_ErrorTopo
global CompareProfilesAreaTopo_Handle
nOut=nargout(CompareProfilesAreaTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesAreaTopo_Handle(varargin{:});
end
