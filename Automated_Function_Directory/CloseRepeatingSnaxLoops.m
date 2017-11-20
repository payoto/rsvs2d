function [varargout]=CloseRepeatingSnaxLoops(varargin)
% include_SnakeSensiv
global CloseRepeatingSnaxLoops_Handle
nOut=nargout(CloseRepeatingSnaxLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CloseRepeatingSnaxLoops_Handle(varargin{:});
end
