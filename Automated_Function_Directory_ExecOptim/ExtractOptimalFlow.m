function [varargout]=ExtractOptimalFlow(varargin)
% OptimisationOutput
global ExtractOptimalFlow_Handle
nOut=nargout(ExtractOptimalFlow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractOptimalFlow_Handle(varargin{:});
end
