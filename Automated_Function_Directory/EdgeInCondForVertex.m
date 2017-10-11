function [varargout]=EdgeInCondForVertex(varargin)
% include_EdgeInformation
global EdgeInCondForVertex_Handle
nOut=nargout(EdgeInCondForVertex_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeInCondForVertex_Handle(varargin{:});
end
