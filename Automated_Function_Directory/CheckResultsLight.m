function [varargout]=CheckResultsLight(varargin)
% include_CheckResultsLight
global CheckResultsLight_Handle
try
nOut=nargout(CheckResultsLight_Handle);
catch
include_CheckResultsLight
nOut=nargout(CheckResultsLight_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckResultsLight_Handle(varargin{:});
end
