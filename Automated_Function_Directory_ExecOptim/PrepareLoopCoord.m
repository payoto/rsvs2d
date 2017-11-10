function [varargout]=PrepareLoopCoord(varargin)
global PrepareLoopCoord_Handle
nOut=nargout(PrepareLoopCoord_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PrepareLoopCoord_Handle(varargin{:});
end
