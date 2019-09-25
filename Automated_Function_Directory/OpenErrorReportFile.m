function [varargout]=OpenErrorReportFile(varargin)
% include_PostProcessing
global OpenErrorReportFile_Handle
try
nOut=nargout(OpenErrorReportFile_Handle);
catch
include_PostProcessing
nOut=nargout(OpenErrorReportFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenErrorReportFile_Handle(varargin{:});
end
