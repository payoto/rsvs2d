function [varargout]=BuildVertexOverflowStruct(varargin)
% include_Optimisation
global BuildVertexOverflowStruct_Handle
nOut=nargout(BuildVertexOverflowStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildVertexOverflowStruct_Handle(varargin{:});
end
