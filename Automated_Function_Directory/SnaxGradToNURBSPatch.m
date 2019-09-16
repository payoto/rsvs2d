function [varargout]=SnaxGradToNURBSPatch(varargin)
% include_NURBSEngine
global SnaxGradToNURBSPatch_Handle
nOut=nargout(SnaxGradToNURBSPatch_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SnaxGradToNURBSPatch_Handle(varargin{:});
end
