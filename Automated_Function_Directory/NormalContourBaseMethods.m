function [varargout]=NormalContourBaseMethods(varargin)
% include_NURBSEngine
global NormalContourBaseMethods_Handle
nOut=nargout(NormalContourBaseMethods_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NormalContourBaseMethods_Handle(varargin{:});
end
