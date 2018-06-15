function [varargout]=PreProcImage(varargin)
% include_Utilities
global PreProcImage_Handle
nOut=nargout(PreProcImage_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PreProcImage_Handle(varargin{:});
end
