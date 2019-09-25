function [varargout]=SplineBasis(varargin)
% include_NURBSEngine
global SplineBasis_Handle
try
nOut=nargout(SplineBasis_Handle);
catch
include_NURBSEngine
nOut=nargout(SplineBasis_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SplineBasis_Handle(varargin{:});
end
