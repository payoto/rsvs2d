function [varargout]=OpenCommentsFile(varargin)
global OpenCommentsFile_Handle
nOut=nargout(OpenCommentsFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenCommentsFile_Handle(varargin{:});
end
