function [varargout]=GoldenSection_func(varargin)
% include_Optimisation
global GoldenSection_func_Handle
try
nOut=nargout(GoldenSection_func_Handle);
catch
include_Optimisation
nOut=nargout(GoldenSection_func_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GoldenSection_func_Handle(varargin{:});
end
