function [varargout]=ExtractSnaxelLoops(varargin)
% include_SnakeParam
global ExtractSnaxelLoops_Handle
nOut=nargout(ExtractSnaxelLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractSnaxelLoops_Handle(varargin{:});
end
