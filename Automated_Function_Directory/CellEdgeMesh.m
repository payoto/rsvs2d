function [varargout]=CellEdgeMesh(varargin)
% include_PostProcessing
global CellEdgeMesh_Handle
nOut=nargout(CellEdgeMesh_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellEdgeMesh_Handle(varargin{:});
end
