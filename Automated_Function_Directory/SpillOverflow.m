function [varargout]=SpillOverflow(varargin)
% include_Optimisation
global SpillOverflow_Handle
try
nOut=nargout(SpillOverflow_Handle);
catch
include_Optimisation
nOut=nargout(SpillOverflow_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SpillOverflow_Handle(varargin{:});
end
