function [varargout]=IndexMod(varargin)
% include_SnakeParam
global IndexMod_Handle
nOut=nargout(IndexMod_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IndexMod_Handle(varargin{:});
end
