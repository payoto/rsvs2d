function [varargout]=OpenBoundaryFile(varargin)
global OpenBoundaryFile_Handle
nOut=nargout(OpenBoundaryFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenBoundaryFile_Handle(varargin{:});
end
