function [varargout]=BuildFillFlow(varargin)
% include_Optimisation
global BuildFillFlow_Handle
try
nOut=nargout(BuildFillFlow_Handle);
catch
include_Optimisation
nOut=nargout(BuildFillFlow_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildFillFlow_Handle(varargin{:});
end
