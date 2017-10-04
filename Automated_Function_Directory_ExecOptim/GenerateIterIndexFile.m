function [varargout]=GenerateIterIndexFile(varargin)
% OptimisationOutput
global GenerateIterIndexFile_Handle
nOut=nargout(GenerateIterIndexFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateIterIndexFile_Handle(varargin{:});
end
