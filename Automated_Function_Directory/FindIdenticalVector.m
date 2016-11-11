function [varargout]=FindIdenticalVector(varargin)
global FindIdenticalVector_Handle
nOut=nargout(FindIdenticalVector_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindIdenticalVector_Handle(varargin{:});
end
