function [varargout]=TrimZerosImage(varargin)
% include_Utilities
global TrimZerosImage_Handle
try
nOut=nargout(TrimZerosImage_Handle);
catch
include_Utilities
nOut=nargout(TrimZerosImage_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimZerosImage_Handle(varargin{:});
end
