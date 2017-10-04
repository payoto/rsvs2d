function [varargout]=OpenOptimumSnakLayFile(varargin)
% OptimisationOutput
global OpenOptimumSnakLayFile_Handle
nOut=nargout(OpenOptimumSnakLayFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenOptimumSnakLayFile_Handle(varargin{:});
end
