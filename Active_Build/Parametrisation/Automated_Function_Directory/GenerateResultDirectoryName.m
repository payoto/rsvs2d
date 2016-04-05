function [varargout]=GenerateResultDirectoryName(varargin)
global GenerateResultDirectoryName_Handle
nOut=nargout(GenerateResultDirectoryName_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateResultDirectoryName_Handle(varargin{:});
end
