function [varargout]=CompressSnaxelChain(varargin)
% include_SnakeParam
global CompressSnaxelChain_Handle
try
nOut=nargout(CompressSnaxelChain_Handle);
catch
include_SnakeParam
nOut=nargout(CompressSnaxelChain_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompressSnaxelChain_Handle(varargin{:});
end
