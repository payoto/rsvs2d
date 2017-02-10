function [varargout]=UpdateStepDir(varargin)
global UpdateStepDir_Handle
nOut=nargout(UpdateStepDir_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=UpdateStepDir_Handle(varargin{:});
end
