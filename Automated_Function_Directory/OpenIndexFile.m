function [varargout]=OpenIndexFile(varargin)
% include_PostProcessing
global OpenIndexFile_Handle
nOut=nargout(OpenIndexFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenIndexFile_Handle(varargin{:});
end
