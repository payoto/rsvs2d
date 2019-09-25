function [varargout]=FindIdenticalVectorOrd(varargin)
% include_Utilities
global FindIdenticalVectorOrd_Handle
try
nOut=nargout(FindIdenticalVectorOrd_Handle);
catch
include_Utilities
nOut=nargout(FindIdenticalVectorOrd_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindIdenticalVectorOrd_Handle(varargin{:});
end
