function [varargout]=CellCentredSnaxelInfo(varargin)
% include_SnakeParam
global CellCentredSnaxelInfo_Handle
nOut=nargout(CellCentredSnaxelInfo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CellCentredSnaxelInfo_Handle(varargin{:});
end
