function [varargout]=GenerateLibCoord(varargin)
global GenerateLibCoord_Handle
nOut=nargout(GenerateLibCoord_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateLibCoord_Handle(varargin{:});
end
