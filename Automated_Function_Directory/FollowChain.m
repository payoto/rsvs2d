function [varargout]=FollowChain(varargin)
% include_EdgeInformation
global FollowChain_Handle
try
nOut=nargout(FollowChain_Handle);
catch
include_EdgeInformation
nOut=nargout(FollowChain_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FollowChain_Handle(varargin{:});
end
