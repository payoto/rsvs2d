function [varargout]=NormalContourBaseMethods(varargin)
% include_NURBSEngine
global NormalContourBaseMethods_Handle
try
nOut=nargout(NormalContourBaseMethods_Handle);
catch
include_NURBSEngine
nOut=nargout(NormalContourBaseMethods_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NormalContourBaseMethods_Handle(varargin{:});
end
