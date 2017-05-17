function [varargout]=GoldenSection_func(varargin)
global GoldenSection_func_Handle
nOut=nargout(GoldenSection_func_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GoldenSection_func_Handle(varargin{:});
end
