function [varargout]=ResizeImage(varargin)
% include_Utilities
global ResizeImage_Handle
try
nOut=nargout(ResizeImage_Handle);
catch
include_Utilities
nOut=nargout(ResizeImage_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ResizeImage_Handle(varargin{:});
end
