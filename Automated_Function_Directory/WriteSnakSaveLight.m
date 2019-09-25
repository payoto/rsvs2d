function [varargout]=WriteSnakSaveLight(varargin)
% include_SnakeParam
global WriteSnakSaveLight_Handle
try
nOut=nargout(WriteSnakSaveLight_Handle);
catch
include_SnakeParam
nOut=nargout(WriteSnakSaveLight_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteSnakSaveLight_Handle(varargin{:});
end
