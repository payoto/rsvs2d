function [varargout]=CheckResults(varargin)
% include_CheckResultsLight
global CheckResults_Handle
try
nOut=nargout(CheckResults_Handle);
catch
include_CheckResultsLight
nOut=nargout(CheckResults_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckResults_Handle(varargin{:});
end
