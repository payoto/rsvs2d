function [varargout]=PrintStart(varargin)
global PrintStart_Handle
nOut=nargout(PrintStart_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PrintStart_Handle(varargin{:});
end
