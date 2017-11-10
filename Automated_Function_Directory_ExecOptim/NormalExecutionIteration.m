function [varargout]=NormalExecutionIteration(varargin)
global NormalExecutionIteration_Handle
nOut=nargout(NormalExecutionIteration_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NormalExecutionIteration_Handle(varargin{:});
end
