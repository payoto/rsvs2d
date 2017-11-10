function [varargout]=CheckIfAnalytical(varargin)
% include_Optimisation
global CheckIfAnalytical_Handle
nOut=nargout(CheckIfAnalytical_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckIfAnalytical_Handle(varargin{:});
end
