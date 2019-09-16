function [varargout]=NormalContourAlternateMethods(varargin)
% include_NURBSEngine
global NormalContourAlternateMethods_Handle
nOut=nargout(NormalContourAlternateMethods_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NormalContourAlternateMethods_Handle(varargin{:});
end
