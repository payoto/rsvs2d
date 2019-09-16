function [varargout]=ExtractContourSnaxel(varargin)
% include_NURBSEngine
global ExtractContourSnaxel_Handle
nOut=nargout(ExtractContourSnaxel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractContourSnaxel_Handle(varargin{:});
end
