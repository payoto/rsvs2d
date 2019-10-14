function [varargout]=OutputAndRunFlowSolve(varargin)
% include_Optimisation
global OutputAndRunFlowSolve_Handle
try
nOut=nargout(OutputAndRunFlowSolve_Handle);
catch
include_Optimisation
nOut=nargout(OutputAndRunFlowSolve_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OutputAndRunFlowSolve_Handle(varargin{:});
end
