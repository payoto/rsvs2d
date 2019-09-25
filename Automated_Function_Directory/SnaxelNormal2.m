function [varargout]=SnaxelNormal2(varargin)
% include_NURBSEngine
global SnaxelNormal2_Handle
try
nOut=nargout(SnaxelNormal2_Handle);
catch
include_NURBSEngine
nOut=nargout(SnaxelNormal2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SnaxelNormal2_Handle(varargin{:});
end
