function [varargout]=CopyFileLinux(varargin)
global CopyFileLinux_Handle
nOut=nargout(CopyFileLinux_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CopyFileLinux_Handle(varargin{:});
end
