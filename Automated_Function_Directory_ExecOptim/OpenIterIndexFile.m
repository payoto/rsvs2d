function [varargout]=OpenIterIndexFile(varargin)
% OptimisationOutput
global OpenIterIndexFile_Handle
nOut=nargout(OpenIterIndexFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenIterIndexFile_Handle(varargin{:});
end
