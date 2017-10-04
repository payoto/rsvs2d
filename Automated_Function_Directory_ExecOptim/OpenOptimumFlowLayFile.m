function [varargout]=OpenOptimumFlowLayFile(varargin)
% OptimisationOutput
global OpenOptimumFlowLayFile_Handle
nOut=nargout(OpenOptimumFlowLayFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenOptimumFlowLayFile_Handle(varargin{:});
end
