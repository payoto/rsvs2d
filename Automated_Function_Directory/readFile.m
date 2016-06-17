function [varargout]=readFile(varargin)
global readFile_Handle
nOut=nargout(readFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=readFile_Handle(varargin{:});
end
