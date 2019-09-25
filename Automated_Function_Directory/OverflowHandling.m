function [varargout]=OverflowHandling(varargin)
% include_Optimisation
global OverflowHandling_Handle
try
nOut=nargout(OverflowHandling_Handle);
catch
include_Optimisation
nOut=nargout(OverflowHandling_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OverflowHandling_Handle(varargin{:});
end
