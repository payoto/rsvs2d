function [varargout]=CompareProfilesAreaSquaredpAreaTopo(varargin)
% InverseDesign_ErrorTopo
global CompareProfilesAreaSquaredpAreaTopo_Handle
nOut=nargout(CompareProfilesAreaSquaredpAreaTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesAreaSquaredpAreaTopo_Handle(varargin{:});
end
