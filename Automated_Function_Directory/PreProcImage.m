function [varargout]=PreProcImage(varargin)
% include_Utilities
global PreProcImage_Handle
try
nOut=nargout(PreProcImage_Handle);
catch
include_Utilities
nOut=nargout(PreProcImage_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PreProcImage_Handle(varargin{:});
end
