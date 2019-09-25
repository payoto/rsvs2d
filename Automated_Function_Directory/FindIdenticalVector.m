function [varargout]=FindIdenticalVector(varargin)
% include_Utilities
global FindIdenticalVector_Handle
try
nOut=nargout(FindIdenticalVector_Handle);
catch
include_Utilities
nOut=nargout(FindIdenticalVector_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindIdenticalVector_Handle(varargin{:});
end
