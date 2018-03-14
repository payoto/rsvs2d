function [varargout]=LoadLoopDataFromMat(varargin)
% InverseDesign_ErrorTopo
global LoadLoopDataFromMat_Handle
nOut=nargout(LoadLoopDataFromMat_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LoadLoopDataFromMat_Handle(varargin{:});
end
