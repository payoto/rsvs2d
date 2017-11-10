function [varargout]=GoldenSection_maxTanDel(varargin)
% include_Optimisation
global GoldenSection_maxTanDel_Handle
nOut=nargout(GoldenSection_maxTanDel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GoldenSection_maxTanDel_Handle(varargin{:});
end
