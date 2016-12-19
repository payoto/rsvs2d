function [varargout]=IterateNoSnake(varargin)
global IterateNoSnake_Handle
nOut=nargout(IterateNoSnake_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IterateNoSnake_Handle(varargin{:});
end
