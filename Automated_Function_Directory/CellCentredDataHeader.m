function [varargout]=CellCentredDataHeader(varargin)
% include_PostProcessing
global CellCentredDataHeader_Handle
try
nOut=nargout(CellCentredDataHeader_Handle);
catch
include_PostProcessing
nOut=nargout(CellCentredDataHeader_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredDataHeader_Handle(varargin{:});
end
