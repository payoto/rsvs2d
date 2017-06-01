function [varargout]=VariableShareZoneString(varargin)
global VariableShareZoneString_Handle
nOut=nargout(VariableShareZoneString_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VariableShareZoneString_Handle(varargin{:});
end
