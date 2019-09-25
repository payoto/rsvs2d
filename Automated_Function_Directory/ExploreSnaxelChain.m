function [varargout]=ExploreSnaxelChain(varargin)
% include_SnakeParam
global ExploreSnaxelChain_Handle
try
nOut=nargout(ExploreSnaxelChain_Handle);
catch
include_SnakeParam
nOut=nargout(ExploreSnaxelChain_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExploreSnaxelChain_Handle(varargin{:});
end
