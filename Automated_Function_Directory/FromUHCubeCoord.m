function [varargout]=FromUHCubeCoord(varargin)
global FromUHCubeCoord_Handle
nOut=nargout(FromUHCubeCoord_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FromUHCubeCoord_Handle(varargin{:});
end
