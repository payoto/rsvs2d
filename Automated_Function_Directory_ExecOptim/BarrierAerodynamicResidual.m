function [varargout]=BarrierAerodynamicResidual(varargin)
global BarrierAerodynamicResidual_Handle
nOut=nargout(BarrierAerodynamicResidual_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BarrierAerodynamicResidual_Handle(varargin{:});
end
