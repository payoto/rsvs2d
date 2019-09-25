function [varargout]=GridRefine_Wrapper(varargin)
% include_Mex_Wrapper
global GridRefine_Wrapper_Handle
try
nOut=nargout(GridRefine_Wrapper_Handle);
catch
include_Mex_Wrapper
nOut=nargout(GridRefine_Wrapper_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GridRefine_Wrapper_Handle(varargin{:});
end
