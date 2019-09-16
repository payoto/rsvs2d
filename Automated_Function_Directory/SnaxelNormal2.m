function [varargout]=SnaxelNormal2(varargin)
% include_NURBSEngine
global SnaxelNormal2_Handle
nOut=nargout(SnaxelNormal2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SnaxelNormal2_Handle(varargin{:});
end
