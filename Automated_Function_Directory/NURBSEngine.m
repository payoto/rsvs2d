function [varargout]=NURBSEngine(varargin)
% include_NURBSEngine
global NURBSEngine_Handle
try
nOut=nargout(NURBSEngine_Handle);
catch
include_NURBSEngine
nOut=nargout(NURBSEngine_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NURBSEngine_Handle(varargin{:});
end
