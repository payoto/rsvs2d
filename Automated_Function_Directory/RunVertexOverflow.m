function [varargout]=RunVertexOverflow(varargin)
% include_Optimisation
global RunVertexOverflow_Handle
try
nOut=nargout(RunVertexOverflow_Handle);
catch
include_Optimisation
nOut=nargout(RunVertexOverflow_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RunVertexOverflow_Handle(varargin{:});
end
