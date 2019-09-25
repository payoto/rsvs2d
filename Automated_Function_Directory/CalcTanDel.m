function [varargout]=CalcTanDel(varargin)
% include_Optimisation
global CalcTanDel_Handle
try
nOut=nargout(CalcTanDel_Handle);
catch
include_Optimisation
nOut=nargout(CalcTanDel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalcTanDel_Handle(varargin{:});
end
