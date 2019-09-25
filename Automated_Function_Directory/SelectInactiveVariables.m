function [varargout]=SelectInactiveVariables(varargin)
% include_Optimisation
global SelectInactiveVariables_Handle
try
nOut=nargout(SelectInactiveVariables_Handle);
catch
include_Optimisation
nOut=nargout(SelectInactiveVariables_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SelectInactiveVariables_Handle(varargin{:});
end
