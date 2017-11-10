function [varargout]=CheckIfGradient(varargin)
% include_Optimisation
global CheckIfGradient_Handle
nOut=nargout(CheckIfGradient_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckIfGradient_Handle(varargin{:});
end
