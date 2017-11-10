function [varargout]=FollowChain(varargin)
% include_EdgeInformation
global FollowChain_Handle
nOut=nargout(FollowChain_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FollowChain_Handle(varargin{:});
end
