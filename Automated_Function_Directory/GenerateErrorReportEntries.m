function [varargout]=GenerateErrorReportEntries(varargin)
% include_PostProcessing
global GenerateErrorReportEntries_Handle
try
nOut=nargout(GenerateErrorReportEntries_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateErrorReportEntries_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateErrorReportEntries_Handle(varargin{:});
end
