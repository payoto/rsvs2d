function [varargout]=StartParallelPool(varargin)
% include_Optimisation
global StartParallelPool_Handle
try
nOut=nargout(StartParallelPool_Handle);
catch
include_Optimisation
nOut=nargout(StartParallelPool_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=StartParallelPool_Handle(varargin{:});
end
