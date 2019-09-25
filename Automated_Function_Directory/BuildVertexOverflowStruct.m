function [varargout]=BuildVertexOverflowStruct(varargin)
% include_Optimisation
global BuildVertexOverflowStruct_Handle
try
nOut=nargout(BuildVertexOverflowStruct_Handle);
catch
include_Optimisation
nOut=nargout(BuildVertexOverflowStruct_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildVertexOverflowStruct_Handle(varargin{:});
end
