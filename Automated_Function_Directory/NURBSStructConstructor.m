function [varargout]=NURBSStructConstructor(varargin)
% include_NURBSEngine
global NURBSStructConstructor_Handle
try
nOut=nargout(NURBSStructConstructor_Handle);
catch
include_NURBSEngine
nOut=nargout(NURBSStructConstructor_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NURBSStructConstructor_Handle(varargin{:});
end
