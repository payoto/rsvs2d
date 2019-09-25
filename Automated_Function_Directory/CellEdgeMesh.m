function [varargout]=CellEdgeMesh(varargin)
% include_PostProcessing
global CellEdgeMesh_Handle
try
nOut=nargout(CellEdgeMesh_Handle);
catch
include_PostProcessing
nOut=nargout(CellEdgeMesh_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellEdgeMesh_Handle(varargin{:});
end
