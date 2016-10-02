function [varargout]=ExecuteOptimisation2(varargin)
global ExecuteOptimisation2_Handle
nOut=nargout(ExecuteOptimisation2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExecuteOptimisation2_Handle(varargin{:});
end
