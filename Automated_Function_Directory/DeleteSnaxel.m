function [varargout]=DeleteSnaxel(varargin)
% include_NURBSEngine
global DeleteSnaxel_Handle
try
nOut=nargout(DeleteSnaxel_Handle);
catch
include_NURBSEngine
nOut=nargout(DeleteSnaxel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DeleteSnaxel_Handle(varargin{:});
end
