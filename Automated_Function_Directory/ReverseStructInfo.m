function [varargout]=ReverseStructInfo(varargin)
% include_Utilities
global ReverseStructInfo_Handle
try
nOut=nargout(ReverseStructInfo_Handle);
catch
include_Utilities
nOut=nargout(ReverseStructInfo_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReverseStructInfo_Handle(varargin{:});
end
