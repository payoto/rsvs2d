function [varargout]=GenerateErrorReportEntries(varargin)
global GenerateErrorReportEntries_Handle
nOut=nargout(GenerateErrorReportEntries_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateErrorReportEntries_Handle(varargin{:});
end
