function [varargout]=WriteToFile(varargin)
% include_PostProcessing
global WriteToFile_Handle
nOut=nargout(WriteToFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteToFile_Handle(varargin{:});
end
