function [varargout]=GoldenSection_maxTanDel(varargin)
% include_Optimisation
global GoldenSection_maxTanDel_Handle
try
nOut=nargout(GoldenSection_maxTanDel_Handle);
catch
include_Optimisation
nOut=nargout(GoldenSection_maxTanDel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GoldenSection_maxTanDel_Handle(varargin{:});
end
