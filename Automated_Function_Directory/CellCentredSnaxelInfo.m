function [varargout]=CellCentredSnaxelInfo(varargin)
% include_SnakeParam
global CellCentredSnaxelInfo_Handle
try
nOut=nargout(CellCentredSnaxelInfo_Handle);
catch
include_SnakeParam
nOut=nargout(CellCentredSnaxelInfo_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredSnaxelInfo_Handle(varargin{:});
end
