function [varargout]=ReadLoopIn(varargin)
% include_Optimisation
global ReadLoopIn_Handle
nOut=nargout(ReadLoopIn_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadLoopIn_Handle(varargin{:});
end
