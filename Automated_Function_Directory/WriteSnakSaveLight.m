function [varargout]=WriteSnakSaveLight(varargin)
% include_SnakeParam
global WriteSnakSaveLight_Handle
nOut=nargout(WriteSnakSaveLight_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteSnakSaveLight_Handle(varargin{:});
end
