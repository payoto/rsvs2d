function [varargout]=SplitSnaxLoops(varargin)
% include_SnakeSensiv
global SplitSnaxLoops_Handle
nOut=nargout(SplitSnaxLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SplitSnaxLoops_Handle(varargin{:});
end
