function [varargout]=NewToOldCell(varargin)
% include_SnakeSensiv
global NewToOldCell_Handle
try
nOut=nargout(NewToOldCell_Handle);
catch
include_SnakeSensiv
nOut=nargout(NewToOldCell_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NewToOldCell_Handle(varargin{:});
end
