function [varargout]=BuildCellConnectivity(varargin)
% include_SnakeSensiv
global BuildCellConnectivity_Handle
nOut=nargout(BuildCellConnectivity_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildCellConnectivity_Handle(varargin{:});
end
