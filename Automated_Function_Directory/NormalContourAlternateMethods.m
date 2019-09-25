function [varargout]=NormalContourAlternateMethods(varargin)
% include_NURBSEngine
global NormalContourAlternateMethods_Handle
try
nOut=nargout(NormalContourAlternateMethods_Handle);
catch
include_NURBSEngine
nOut=nargout(NormalContourAlternateMethods_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NormalContourAlternateMethods_Handle(varargin{:});
end
