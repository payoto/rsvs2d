function [varargout]=BuildNewRestartFrac(varargin)
global BuildNewRestartFrac_Handle
nOut=nargout(BuildNewRestartFrac_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildNewRestartFrac_Handle(varargin{:});
end
