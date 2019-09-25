function [varargout]=CloseRepeatingSnaxLoops(varargin)
% include_SnakeSensiv
global CloseRepeatingSnaxLoops_Handle
try
nOut=nargout(CloseRepeatingSnaxLoops_Handle);
catch
include_SnakeSensiv
nOut=nargout(CloseRepeatingSnaxLoops_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CloseRepeatingSnaxLoops_Handle(varargin{:});
end
