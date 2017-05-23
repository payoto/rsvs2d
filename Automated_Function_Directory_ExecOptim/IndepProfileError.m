function [varargout]=IndepProfileError(varargin)
global IndepProfileError_Handle
nOut=nargout(IndepProfileError_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IndepProfileError_Handle(varargin{:});
end
