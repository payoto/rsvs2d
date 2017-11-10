function [varargout]=PrintEnd(varargin)
global PrintEnd_Handle
nOut=nargout(PrintEnd_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PrintEnd_Handle(varargin{:});
end
