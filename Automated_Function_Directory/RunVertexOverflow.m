function [varargout]=RunVertexOverflow(varargin)
global RunVertexOverflow_Handle
nOut=nargout(RunVertexOverflow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RunVertexOverflow_Handle(varargin{:});
end
