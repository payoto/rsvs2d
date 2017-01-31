function [varargout]=GenerateSnakStartLoop(varargin)
global GenerateSnakStartLoop_Handle
nOut=nargout(GenerateSnakStartLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateSnakStartLoop_Handle(varargin{:});
end
