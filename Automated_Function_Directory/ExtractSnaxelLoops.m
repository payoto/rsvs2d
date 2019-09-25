function [varargout]=ExtractSnaxelLoops(varargin)
% include_SnakeParam
global ExtractSnaxelLoops_Handle
try
nOut=nargout(ExtractSnaxelLoops_Handle);
catch
include_SnakeParam
nOut=nargout(ExtractSnaxelLoops_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractSnaxelLoops_Handle(varargin{:});
end
