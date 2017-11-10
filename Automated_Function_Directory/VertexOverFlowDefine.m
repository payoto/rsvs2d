function [varargout]=VertexOverFlowDefine(varargin)
% include_Optimisation
global VertexOverFlowDefine_Handle
nOut=nargout(VertexOverFlowDefine_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VertexOverFlowDefine_Handle(varargin{:});
end
