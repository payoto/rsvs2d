function [varargout]=RecursiveStringGeneration(varargin)
% include_PostProcessing
global RecursiveStringGeneration_Handle
nOut=nargout(RecursiveStringGeneration_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RecursiveStringGeneration_Handle(varargin{:});
end
