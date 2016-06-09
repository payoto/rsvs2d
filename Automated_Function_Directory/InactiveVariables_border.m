function [varargout]=InactiveVariables_border(varargin)
global InactiveVariables_border_Handle
nOut=nargout(InactiveVariables_border_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InactiveVariables_border_Handle(varargin{:});
end
