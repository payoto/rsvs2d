function [varargout]=TrimZeros(varargin)
% include_Utilities
global TrimZeros_Handle
try
nOut=nargout(TrimZeros_Handle);
catch
include_Utilities
nOut=nargout(TrimZeros_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimZeros_Handle(varargin{:});
end
