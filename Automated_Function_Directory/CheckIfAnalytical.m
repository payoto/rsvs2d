function [varargout]=CheckIfAnalytical(varargin)
% include_Optimisation
global CheckIfAnalytical_Handle
try
nOut=nargout(CheckIfAnalytical_Handle);
catch
include_Optimisation
nOut=nargout(CheckIfAnalytical_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckIfAnalytical_Handle(varargin{:});
end
