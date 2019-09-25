function [varargout]=ComputeArcArea(varargin)
% include_NURBSEngine
global ComputeArcArea_Handle
try
nOut=nargout(ComputeArcArea_Handle);
catch
include_NURBSEngine
nOut=nargout(ComputeArcArea_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeArcArea_Handle(varargin{:});
end
