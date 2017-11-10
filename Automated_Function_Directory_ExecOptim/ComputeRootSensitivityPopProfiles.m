function [varargout]=ComputeRootSensitivityPopProfiles(varargin)
global ComputeRootSensitivityPopProfiles_Handle
nOut=nargout(ComputeRootSensitivityPopProfiles_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeRootSensitivityPopProfiles_Handle(varargin{:});
end
