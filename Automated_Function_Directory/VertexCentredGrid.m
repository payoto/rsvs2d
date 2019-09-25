function [varargout]=VertexCentredGrid(varargin)
% include_SnakeParam
global VertexCentredGrid_Handle
try
nOut=nargout(VertexCentredGrid_Handle);
catch
include_SnakeParam
nOut=nargout(VertexCentredGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VertexCentredGrid_Handle(varargin{:});
end
