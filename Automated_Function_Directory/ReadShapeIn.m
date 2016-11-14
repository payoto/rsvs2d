function [varargout]=ReadShapeIn(varargin)
global ReadShapeIn_Handle
nOut=nargout(ReadShapeIn_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadShapeIn_Handle(varargin{:});
end
