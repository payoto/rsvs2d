function [varargout]=ExtractUpperLowerSurface(varargin)
% include_Optimisation
global ExtractUpperLowerSurface_Handle
try
nOut=nargout(ExtractUpperLowerSurface_Handle);
catch
include_Optimisation
nOut=nargout(ExtractUpperLowerSurface_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractUpperLowerSurface_Handle(varargin{:});
end
