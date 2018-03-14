function [varargout]=CompareProfilesNormDistTopo(varargin)
% InverseDesign_ErrorTopo
global CompareProfilesNormDistTopo_Handle
nOut=nargout(CompareProfilesNormDistTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesNormDistTopo_Handle(varargin{:});
end
