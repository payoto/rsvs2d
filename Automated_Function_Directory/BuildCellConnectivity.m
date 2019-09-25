function [varargout]=BuildCellConnectivity(varargin)
% include_SnakeSensiv
global BuildCellConnectivity_Handle
try
nOut=nargout(BuildCellConnectivity_Handle);
catch
include_SnakeSensiv
nOut=nargout(BuildCellConnectivity_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildCellConnectivity_Handle(varargin{:});
end
