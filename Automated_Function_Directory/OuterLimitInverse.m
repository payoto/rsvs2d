function [varargout]=OuterLimitInverse(varargin)
% include_Optimisation
global OuterLimitInverse_Handle
try
nOut=nargout(OuterLimitInverse_Handle);
catch
include_Optimisation
nOut=nargout(OuterLimitInverse_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OuterLimitInverse_Handle(varargin{:});
end
