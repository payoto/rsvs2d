function [varargout]=VertexOverFlowDefine(varargin)
% include_Optimisation
global VertexOverFlowDefine_Handle
try
nOut=nargout(VertexOverFlowDefine_Handle);
catch
include_Optimisation
nOut=nargout(VertexOverFlowDefine_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VertexOverFlowDefine_Handle(varargin{:});
end
