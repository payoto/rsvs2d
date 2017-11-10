function [varargout]=OpenOptimLayFile(varargin)
% OptimisationOutput
global OpenOptimLayFile_Handle
nOut=nargout(OpenOptimLayFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenOptimLayFile_Handle(varargin{:});
end
