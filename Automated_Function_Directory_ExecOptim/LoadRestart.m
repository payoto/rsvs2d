function [varargout]=LoadRestart(varargin)
global LoadRestart_Handle
nOut=nargout(LoadRestart_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LoadRestart_Handle(varargin{:});
end
