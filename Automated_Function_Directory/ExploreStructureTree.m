function [varargout]=ExploreStructureTree(varargin)
% include_Utilities
global ExploreStructureTree_Handle
try
nOut=nargout(ExploreStructureTree_Handle);
catch
include_Utilities
nOut=nargout(ExploreStructureTree_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExploreStructureTree_Handle(varargin{:});
end
