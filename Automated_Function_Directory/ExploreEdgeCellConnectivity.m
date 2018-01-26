function [varargout]=ExploreEdgeCellConnectivity(varargin)
% include_SnakeParam
global ExploreEdgeCellConnectivity_Handle
nOut=nargout(ExploreEdgeCellConnectivity_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExploreEdgeCellConnectivity_Handle(varargin{:});
end
