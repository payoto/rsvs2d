function [varargout]=CellCentredGrid(varargin)
global CellCentredGrid_Handle
nOut=nargout(CellCentredGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredGrid_Handle(varargin{:});
end
