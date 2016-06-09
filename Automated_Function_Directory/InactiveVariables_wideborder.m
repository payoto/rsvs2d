function [varargout]=InactiveVariables_wideborder(varargin)
global InactiveVariables_wideborder_Handle
nOut=nargout(InactiveVariables_wideborder_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InactiveVariables_wideborder_Handle(varargin{:});
end
