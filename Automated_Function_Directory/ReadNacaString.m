function [varargout]=ReadNacaString(varargin)
global ReadNacaString_Handle
nOut=nargout(ReadNacaString_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadNacaString_Handle(varargin{:});
end
