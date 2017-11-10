function [varargout]=IdentifiedCrossedEdge(varargin)
global IdentifiedCrossedEdge_Handle
nOut=nargout(IdentifiedCrossedEdge_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IdentifiedCrossedEdge_Handle(varargin{:});
end
