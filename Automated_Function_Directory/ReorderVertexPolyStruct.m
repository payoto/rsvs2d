function [varargout]=ReorderVertexPolyStruct(varargin)
% include_PostProcessing
global ReorderVertexPolyStruct_Handle
nOut=nargout(ReorderVertexPolyStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReorderVertexPolyStruct_Handle(varargin{:});
end
