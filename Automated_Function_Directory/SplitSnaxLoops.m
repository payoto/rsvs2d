function [varargout]=SplitSnaxLoops(varargin)
% include_SnakeSensiv
global SplitSnaxLoops_Handle
try
nOut=nargout(SplitSnaxLoops_Handle);
catch
include_SnakeSensiv
nOut=nargout(SplitSnaxLoops_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SplitSnaxLoops_Handle(varargin{:});
end
