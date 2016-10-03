function [varargout]=FindDir(varargin)
global FindDir_Handle
nOut=nargout(FindDir_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindDir_Handle(varargin{:});
end
