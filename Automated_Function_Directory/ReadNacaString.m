function [varargout]=ReadNacaString(varargin)
% include_Optimisation
global ReadNacaString_Handle
try
nOut=nargout(ReadNacaString_Handle);
catch
include_Optimisation
nOut=nargout(ReadNacaString_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadNacaString_Handle(varargin{:});
end
