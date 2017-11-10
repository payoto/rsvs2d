function [varargout]=UpdateSupportOptim(varargin)
global UpdateSupportOptim_Handle
nOut=nargout(UpdateSupportOptim_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=UpdateSupportOptim_Handle(varargin{:});
end
