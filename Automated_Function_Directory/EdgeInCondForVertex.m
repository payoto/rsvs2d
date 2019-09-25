function [varargout]=EdgeInCondForVertex(varargin)
% include_EdgeInformation
global EdgeInCondForVertex_Handle
try
nOut=nargout(EdgeInCondForVertex_Handle);
catch
include_EdgeInformation
nOut=nargout(EdgeInCondForVertex_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeInCondForVertex_Handle(varargin{:});
end
