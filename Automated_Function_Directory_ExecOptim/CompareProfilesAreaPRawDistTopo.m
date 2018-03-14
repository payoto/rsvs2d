function [varargout]=CompareProfilesAreaPRawDistTopo(varargin)
% InverseDesign_ErrorTopo
global CompareProfilesAreaPRawDistTopo_Handle
nOut=nargout(CompareProfilesAreaPRawDistTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesAreaPRawDistTopo_Handle(varargin{:});
end
