function [varargout]=CalcTanDel(varargin)
global CalcTanDel_Handle
nOut=nargout(CalcTanDel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalcTanDel_Handle(varargin{:});
end
