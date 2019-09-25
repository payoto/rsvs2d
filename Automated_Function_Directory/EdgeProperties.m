function [varargout]=EdgeProperties(varargin)
% include_EdgeInformation
global EdgeProperties_Handle
try
nOut=nargout(EdgeProperties_Handle);
catch
include_EdgeInformation
nOut=nargout(EdgeProperties_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeProperties_Handle(varargin{:});
end
