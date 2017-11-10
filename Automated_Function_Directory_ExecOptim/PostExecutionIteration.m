function [varargout]=PostExecutionIteration(varargin)
global PostExecutionIteration_Handle
nOut=nargout(PostExecutionIteration_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PostExecutionIteration_Handle(varargin{:});
end
