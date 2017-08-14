function [varargout]=VariableShareZoneString(varargin)
% include_PostProcessing
global VariableShareZoneString_Handle
nOut=nargout(VariableShareZoneString_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VariableShareZoneString_Handle(varargin{:});
end
