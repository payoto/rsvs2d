function [varargout]=ReorderVertexPolyStruct(varargin)
% include_PostProcessing
global ReorderVertexPolyStruct_Handle
try
nOut=nargout(ReorderVertexPolyStruct_Handle);
catch
include_PostProcessing
nOut=nargout(ReorderVertexPolyStruct_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReorderVertexPolyStruct_Handle(varargin{:});
end
