function [varargout]=DataToString(varargin)
global DataToString_Handle
nOut=nargout(DataToString_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DataToString_Handle(varargin{:});
end
