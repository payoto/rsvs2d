function [varargout]=SnaxGradToNURBSPatch(varargin)
% include_NURBSEngine
global SnaxGradToNURBSPatch_Handle
try
nOut=nargout(SnaxGradToNURBSPatch_Handle);
catch
include_NURBSEngine
nOut=nargout(SnaxGradToNURBSPatch_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SnaxGradToNURBSPatch_Handle(varargin{:});
end
