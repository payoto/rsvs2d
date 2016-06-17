function [varargout]=GenerateErrorReportFile(varargin)
global GenerateErrorReportFile_Handle
nOut=nargout(GenerateErrorReportFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateErrorReportFile_Handle(varargin{:});
end
