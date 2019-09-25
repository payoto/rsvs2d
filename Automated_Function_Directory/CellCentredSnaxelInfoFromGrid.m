function [varargout]=CellCentredSnaxelInfoFromGrid(varargin)
% include_NURBSEngine
global CellCentredSnaxelInfoFromGrid_Handle
try
nOut=nargout(CellCentredSnaxelInfoFromGrid_Handle);
catch
include_NURBSEngine
nOut=nargout(CellCentredSnaxelInfoFromGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredSnaxelInfoFromGrid_Handle(varargin{:});
end
