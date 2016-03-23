function [varargout]=CopyDiary(varargin)
global CopyDiary_Handle
nOut=nargout(CopyDiary_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CopyDiary_Handle(varargin{:});
end
