function [varargout]=CellCentredSnaxelInfoFromGrid(varargin)
% include_NURBSEngine
global CellCentredSnaxelInfoFromGrid_Handle
nOut=nargout(CellCentredSnaxelInfoFromGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredSnaxelInfoFromGrid_Handle(varargin{:});
end
