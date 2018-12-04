function [varargout]=ConcatenateTecplotFile(varargin)
% include_PostProcessing
global ConcatenateTecplotFile_Handle
nOut=nargout(ConcatenateTecplotFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConcatenateTecplotFile_Handle(varargin{:});
end
