function [varargout]=EditVariablePLT_FELINESEG(varargin)
% OptimisationOutput
global EditVariablePLT_FELINESEG_Handle
nOut=nargout(EditVariablePLT_FELINESEG_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EditVariablePLT_FELINESEG_Handle(varargin{:});
end
