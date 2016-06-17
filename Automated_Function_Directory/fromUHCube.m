function [varargout]=fromUHCube(varargin)
global fromUHCube_Handle
nOut=nargout(fromUHCube_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=fromUHCube_Handle(varargin{:});
end
