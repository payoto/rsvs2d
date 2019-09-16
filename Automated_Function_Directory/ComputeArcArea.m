function [varargout]=ComputeArcArea(varargin)
% include_NURBSEngine
global ComputeArcArea_Handle
nOut=nargout(ComputeArcArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeArcArea_Handle(varargin{:});
end
