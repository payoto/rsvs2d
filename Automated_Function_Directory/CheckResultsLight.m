function [varargout]=CheckResultsLight(varargin)
% include_CheckResultsLight
global CheckResultsLight_Handle
nOut=nargout(CheckResultsLight_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckResultsLight_Handle(varargin{:});
end
