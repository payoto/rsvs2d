function [varargout]=GenerateSummaryLine(varargin)
% include_Validation
global GenerateSummaryLine_Handle
try
nOut=nargout(GenerateSummaryLine_Handle);
catch
include_Validation
nOut=nargout(GenerateSummaryLine_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateSummaryLine_Handle(varargin{:});
end
