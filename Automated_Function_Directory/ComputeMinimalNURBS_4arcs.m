function [varargout]=ComputeMinimalNURBS_4arcs(varargin)
% include_NURBSEngine
global ComputeMinimalNURBS_4arcs_Handle
nOut=nargout(ComputeMinimalNURBS_4arcs_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeMinimalNURBS_4arcs_Handle(varargin{:});
end
