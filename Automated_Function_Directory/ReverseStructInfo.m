function [varargout]=ReverseStructInfo(varargin)
% include_Utilities
global ReverseStructInfo_Handle
nOut=nargout(ReverseStructInfo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReverseStructInfo_Handle(varargin{:});
end
