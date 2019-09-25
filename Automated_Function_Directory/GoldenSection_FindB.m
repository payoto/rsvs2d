function [varargout]=GoldenSection_FindB(varargin)
% include_Optimisation
global GoldenSection_FindB_Handle
try
nOut=nargout(GoldenSection_FindB_Handle);
catch
include_Optimisation
nOut=nargout(GoldenSection_FindB_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GoldenSection_FindB_Handle(varargin{:});
end
