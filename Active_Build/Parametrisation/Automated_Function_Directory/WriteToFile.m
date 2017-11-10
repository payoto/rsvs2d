function [varargout]=WriteToFile(varargin)
global WriteToFile_Handle
nOut=nargout(WriteToFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteToFile_Handle(varargin{:});
end
