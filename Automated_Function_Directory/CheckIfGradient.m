function [varargout]=CheckIfGradient(varargin)
% include_Optimisation
global CheckIfGradient_Handle
try
nOut=nargout(CheckIfGradient_Handle);
catch
include_Optimisation
nOut=nargout(CheckIfGradient_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckIfGradient_Handle(varargin{:});
end
