function [varargout]=GridInit_Wrapper(varargin)
global GridInit_Wrapper_Handle
nOut=nargout(GridInit_Wrapper_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GridInit_Wrapper_Handle(varargin{:});
end
