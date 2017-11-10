function [varargout]=GridRefine_Wrapper(varargin)
global GridRefine_Wrapper_Handle
nOut=nargout(GridRefine_Wrapper_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GridRefine_Wrapper_Handle(varargin{:});
end
