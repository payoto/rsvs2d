function [varargout]=GenerateSummaryLine(varargin)
% include_Validation
global GenerateSummaryLine_Handle
nOut=nargout(GenerateSummaryLine_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateSummaryLine_Handle(varargin{:});
end
