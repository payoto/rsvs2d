function [varargout]=InitOptimalFlowOutput(varargin)
% OptimisationOutput
global InitOptimalFlowOutput_Handle
nOut=nargout(InitOptimalFlowOutput_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitOptimalFlowOutput_Handle(varargin{:});
end
