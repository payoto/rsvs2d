function [varargout]=fromUHCubeStruct(varargin)
global fromUHCubeStruct_Handle
nOut=nargout(fromUHCubeStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=fromUHCubeStruct_Handle(varargin{:});
end
