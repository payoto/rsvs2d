function [varargout]=Loop2CellEdgeMesh(varargin)
% include_PostProcessing
global Loop2CellEdgeMesh_Handle
try
nOut=nargout(Loop2CellEdgeMesh_Handle);
catch
include_PostProcessing
nOut=nargout(Loop2CellEdgeMesh_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=Loop2CellEdgeMesh_Handle(varargin{:});
end
