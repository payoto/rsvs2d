function [varargout]=NewToOldCell(varargin)
% include_SnakeSensiv
global NewToOldCell_Handle
nOut=nargout(NewToOldCell_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NewToOldCell_Handle(varargin{:});
end
