function [varargout]=OpenErrorReportFile(varargin)
% include_PostProcessing
global OpenErrorReportFile_Handle
nOut=nargout(OpenErrorReportFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenErrorReportFile_Handle(varargin{:});
end
