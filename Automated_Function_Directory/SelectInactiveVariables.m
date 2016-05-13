function [varargout]=SelectInactiveVariables(varargin)
global SelectInactiveVariables_Handle
nOut=nargout(SelectInactiveVariables_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SelectInactiveVariables_Handle(varargin{:});
end
