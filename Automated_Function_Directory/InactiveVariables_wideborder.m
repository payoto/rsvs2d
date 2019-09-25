function [varargout]=InactiveVariables_wideborder(varargin)
% include_Optimisation
global InactiveVariables_wideborder_Handle
try
nOut=nargout(InactiveVariables_wideborder_Handle);
catch
include_Optimisation
nOut=nargout(InactiveVariables_wideborder_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InactiveVariables_wideborder_Handle(varargin{:});
end
