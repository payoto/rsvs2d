function [varargout]=CompareProfilesAreaRawDistTopo(varargin)
% InverseDesign_ErrorTopo
global CompareProfilesAreaRawDistTopo_Handle
nOut=nargout(CompareProfilesAreaRawDistTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesAreaRawDistTopo_Handle(varargin{:});
end
