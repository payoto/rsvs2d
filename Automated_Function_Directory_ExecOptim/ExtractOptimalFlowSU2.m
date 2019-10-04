function [varargout]=ExtractOptimalFlowSU2(varargin)
% OptimisationOutput
global ExtractOptimalFlowSU2_Handle
try
nOut=nargout(ExtractOptimalFlowSU2_Handle);
catch
OptimisationOutput
nOut=nargout(ExtractOptimalFlowSU2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractOptimalFlowSU2_Handle(varargin{:});
end
