function [varargout]=ReadShapeIn(varargin)
% include_Optimisation
global ReadShapeIn_Handle
try
nOut=nargout(ReadShapeIn_Handle);
catch
include_Optimisation
nOut=nargout(ReadShapeIn_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadShapeIn_Handle(varargin{:});
end
