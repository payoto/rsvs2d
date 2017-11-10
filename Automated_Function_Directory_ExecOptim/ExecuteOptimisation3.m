function [varargout]=ExecuteOptimisation3(varargin)
global ExecuteOptimisation3_Handle
nOut=nargout(ExecuteOptimisation3_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExecuteOptimisation3_Handle(varargin{:});
end
