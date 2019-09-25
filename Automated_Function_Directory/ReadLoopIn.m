function [varargout]=ReadLoopIn(varargin)
% include_Optimisation
global ReadLoopIn_Handle
try
nOut=nargout(ReadLoopIn_Handle);
catch
include_Optimisation
nOut=nargout(ReadLoopIn_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadLoopIn_Handle(varargin{:});
end
