function [varargout]=OpenBoundaryFile(varargin)
% include_PostProcessing
global OpenBoundaryFile_Handle
try
nOut=nargout(OpenBoundaryFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenBoundaryFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenBoundaryFile_Handle(varargin{:});
end
