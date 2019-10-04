function [varargout]=ExtractOptimalFlowSU2Multi(varargin)
% OptimisationOutput
global ExtractOptimalFlowSU2Multi_Handle
try
nOut=nargout(ExtractOptimalFlowSU2Multi_Handle);
catch
OptimisationOutput
nOut=nargout(ExtractOptimalFlowSU2Multi_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractOptimalFlowSU2Multi_Handle(varargin{:});
end
