function [varargout]=ReadShapeIn(varargin)
% include_Optimisation
global ReadShapeIn_Handle
nOut=nargout(ReadShapeIn_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadShapeIn_Handle(varargin{:});
end
