function [varargout]=RescaleDesVarNoFill(varargin)
% include_Optimisation
global RescaleDesVarNoFill_Handle
try
nOut=nargout(RescaleDesVarNoFill_Handle);
catch
include_Optimisation
nOut=nargout(RescaleDesVarNoFill_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RescaleDesVarNoFill_Handle(varargin{:});
end
