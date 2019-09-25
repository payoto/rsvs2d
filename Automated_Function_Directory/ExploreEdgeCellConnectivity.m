function [varargout]=ExploreEdgeCellConnectivity(varargin)
% include_SnakeParam
global ExploreEdgeCellConnectivity_Handle
try
nOut=nargout(ExploreEdgeCellConnectivity_Handle);
catch
include_SnakeParam
nOut=nargout(ExploreEdgeCellConnectivity_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExploreEdgeCellConnectivity_Handle(varargin{:});
end
