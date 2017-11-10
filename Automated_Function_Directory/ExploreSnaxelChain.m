function [varargout]=ExploreSnaxelChain(varargin)
% include_SnakeParam
global ExploreSnaxelChain_Handle
nOut=nargout(ExploreSnaxelChain_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExploreSnaxelChain_Handle(varargin{:});
end
