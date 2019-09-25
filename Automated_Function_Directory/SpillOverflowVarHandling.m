function [varargout]=SpillOverflowVarHandling(varargin)
% include_Optimisation
global SpillOverflowVarHandling_Handle
try
nOut=nargout(SpillOverflowVarHandling_Handle);
catch
include_Optimisation
nOut=nargout(SpillOverflowVarHandling_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SpillOverflowVarHandling_Handle(varargin{:});
end
