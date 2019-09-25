function [varargout]=GenerateErrorReportFile(varargin)
% include_PostProcessing
global GenerateErrorReportFile_Handle
try
nOut=nargout(GenerateErrorReportFile_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateErrorReportFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateErrorReportFile_Handle(varargin{:});
end
