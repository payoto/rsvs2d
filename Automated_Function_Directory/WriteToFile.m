function [varargout]=WriteToFile(varargin)
% include_PostProcessing
global WriteToFile_Handle
try
nOut=nargout(WriteToFile_Handle);
catch
include_PostProcessing
nOut=nargout(WriteToFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteToFile_Handle(varargin{:});
end
