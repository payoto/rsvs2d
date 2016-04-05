function [varargout]=MakeCommentsFile(varargin)
global MakeCommentsFile_Handle
nOut=nargout(MakeCommentsFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCommentsFile_Handle(varargin{:});
end
