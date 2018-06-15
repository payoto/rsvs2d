function [varargout]=ResizeImage(varargin)
% include_Utilities
global ResizeImage_Handle
nOut=nargout(ResizeImage_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ResizeImage_Handle(varargin{:});
end
