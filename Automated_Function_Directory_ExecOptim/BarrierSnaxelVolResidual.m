function [varargout]=BarrierSnaxelVolResidual(varargin)
global BarrierSnaxelVolResidual_Handle
nOut=nargout(BarrierSnaxelVolResidual_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BarrierSnaxelVolResidual_Handle(varargin{:});
end
