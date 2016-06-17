function [varargout]=ToUHCubeCoord(varargin)
global ToUHCubeCoord_Handle
nOut=nargout(ToUHCubeCoord_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ToUHCubeCoord_Handle(varargin{:});
end
