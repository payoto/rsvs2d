function [varargout]=BuildNurbsLoops(varargin)
% include_NURBSEngine
global BuildNurbsLoops_Handle
nOut=nargout(BuildNurbsLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildNurbsLoops_Handle(varargin{:});
end
