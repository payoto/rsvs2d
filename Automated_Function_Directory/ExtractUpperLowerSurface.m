function [varargout]=ExtractUpperLowerSurface(varargin)
% include_Optimisation
global ExtractUpperLowerSurface_Handle
nOut=nargout(ExtractUpperLowerSurface_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractUpperLowerSurface_Handle(varargin{:});
end
