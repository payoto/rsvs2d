function [varargout]=CellCentredDataHeader(varargin)
global CellCentredDataHeader_Handle
nOut=nargout(CellCentredDataHeader_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredDataHeader_Handle(varargin{:});
end
