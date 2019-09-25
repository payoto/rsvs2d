function [varargout]=SetVariables(varargin)
% include_Utilities
global SetVariables_Handle
try
nOut=nargout(SetVariables_Handle);
catch
include_Utilities
nOut=nargout(SetVariables_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SetVariables_Handle(varargin{:});
end
