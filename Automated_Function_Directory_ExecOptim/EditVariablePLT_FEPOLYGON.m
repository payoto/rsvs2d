function [varargout]=EditVariablePLT_FEPOLYGON(varargin)
% OptimisationOutput
global EditVariablePLT_FEPOLYGON_Handle
nOut=nargout(EditVariablePLT_FEPOLYGON_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EditVariablePLT_FEPOLYGON_Handle(varargin{:});
end
