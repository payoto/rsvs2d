function [varargout]=VertexCentredGrid(varargin)
% include_SnakeParam
global VertexCentredGrid_Handle
nOut=nargout(VertexCentredGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VertexCentredGrid_Handle(varargin{:});
end
