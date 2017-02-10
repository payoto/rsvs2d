function [varargout]=OptimisationDebug(varargin)
global OptimisationDebug_Handle
nOut=nargout(OptimisationDebug_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationDebug_Handle(varargin{:});
end
