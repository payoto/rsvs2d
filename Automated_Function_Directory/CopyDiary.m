function [varargout]=CopyDiary(varargin)
% include_PostProcessing
global CopyDiary_Handle
try
nOut=nargout(CopyDiary_Handle);
catch
include_PostProcessing
nOut=nargout(CopyDiary_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CopyDiary_Handle(varargin{:});
end
