function [varargout]=InitialiseAeroshell(varargin)
global InitialiseAeroshell_Handle
nOut=nargout(InitialiseAeroshell_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseAeroshell_Handle(varargin{:});
end
