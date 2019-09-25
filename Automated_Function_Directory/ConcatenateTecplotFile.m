function [varargout]=ConcatenateTecplotFile(varargin)
% include_PostProcessing
global ConcatenateTecplotFile_Handle
try
nOut=nargout(ConcatenateTecplotFile_Handle);
catch
include_PostProcessing
nOut=nargout(ConcatenateTecplotFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConcatenateTecplotFile_Handle(varargin{:});
end
