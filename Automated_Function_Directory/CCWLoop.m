function [varargout]=CCWLoop(varargin)
% include_SnakeParam
global CCWLoop_Handle
try
nOut=nargout(CCWLoop_Handle);
catch
include_SnakeParam
nOut=nargout(CCWLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CCWLoop_Handle(varargin{:});
end
