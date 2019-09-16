function [varargout]=DeleteSnaxel(varargin)
% include_NURBSEngine
global DeleteSnaxel_Handle
nOut=nargout(DeleteSnaxel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DeleteSnaxel_Handle(varargin{:});
end
