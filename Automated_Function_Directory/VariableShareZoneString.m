function [varargout]=VariableShareZoneString(varargin)
% include_PostProcessing
global VariableShareZoneString_Handle
try
nOut=nargout(VariableShareZoneString_Handle);
catch
include_PostProcessing
nOut=nargout(VariableShareZoneString_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VariableShareZoneString_Handle(varargin{:});
end
