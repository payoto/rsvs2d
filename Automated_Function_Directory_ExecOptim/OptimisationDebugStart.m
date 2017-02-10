function [varargout]=OptimisationDebugStart(varargin)
global OptimisationDebugStart_Handle
nOut=nargout(OptimisationDebugStart_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationDebugStart_Handle(varargin{:});
end
