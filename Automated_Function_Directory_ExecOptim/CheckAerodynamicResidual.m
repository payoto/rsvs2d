function [varargout]=CheckAerodynamicResidual(varargin)
global CheckAerodynamicResidual_Handle
nOut=nargout(CheckAerodynamicResidual_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckAerodynamicResidual_Handle(varargin{:});
end
