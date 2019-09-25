function [varargout]=IndexMod(varargin)
% include_SnakeParam
global IndexMod_Handle
try
nOut=nargout(IndexMod_Handle);
catch
include_SnakeParam
nOut=nargout(IndexMod_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IndexMod_Handle(varargin{:});
end
