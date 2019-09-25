function [varargout]=ExtractContourSnaxel(varargin)
% include_NURBSEngine
global ExtractContourSnaxel_Handle
try
nOut=nargout(ExtractContourSnaxel_Handle);
catch
include_NURBSEngine
nOut=nargout(ExtractContourSnaxel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractContourSnaxel_Handle(varargin{:});
end
