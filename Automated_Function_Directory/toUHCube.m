function [varargout]=toUHCube(varargin)
global toUHCube_Handle
nOut=nargout(toUHCube_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=toUHCube_Handle(varargin{:});
end
