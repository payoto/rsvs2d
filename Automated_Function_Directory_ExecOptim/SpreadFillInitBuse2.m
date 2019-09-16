function [varargout]=SpreadFillInitBuse2(varargin)
% ExecuteOptimisation
global SpreadFillInitBuse2_Handle
nOut=nargout(SpreadFillInitBuse2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SpreadFillInitBuse2_Handle(varargin{:});
end
