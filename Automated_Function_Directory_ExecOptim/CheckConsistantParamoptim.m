function [varargout]=CheckConsistantParamoptim(varargin)
% ExecuteOptimisation
global CheckConsistantParamoptim_Handle
nOut=nargout(CheckConsistantParamoptim_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckConsistantParamoptim_Handle(varargin{:});
end
