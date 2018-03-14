function [varargout]=LoadLoopData(varargin)
% InverseDesign_ErrorTopo
global LoadLoopData_Handle
nOut=nargout(LoadLoopData_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LoadLoopData_Handle(varargin{:});
end
