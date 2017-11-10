function [varargout]=ExploreStructureTree(varargin)
global ExploreStructureTree_Handle
nOut=nargout(ExploreStructureTree_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExploreStructureTree_Handle(varargin{:});
end
