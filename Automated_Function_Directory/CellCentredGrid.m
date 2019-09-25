function [varargout]=CellCentredGrid(varargin)
% include_SnakeParam
global CellCentredGrid_Handle
try
nOut=nargout(CellCentredGrid_Handle);
catch
include_SnakeParam
nOut=nargout(CellCentredGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredGrid_Handle(varargin{:});
end
