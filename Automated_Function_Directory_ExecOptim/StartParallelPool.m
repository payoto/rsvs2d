function [varargout]=StartParallelPool(varargin)
global StartParallelPool_Handle
nOut=nargout(StartParallelPool_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=StartParallelPool_Handle(varargin{:});
end
