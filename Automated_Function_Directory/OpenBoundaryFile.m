function [varargout]=OpenBoundaryFile(varargin)
% include_PostProcessing
global OpenBoundaryFile_Handle
nOut=nargout(OpenBoundaryFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenBoundaryFile_Handle(varargin{:});
end
