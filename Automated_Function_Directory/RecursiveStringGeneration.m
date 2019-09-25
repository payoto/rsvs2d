function [varargout]=RecursiveStringGeneration(varargin)
% include_PostProcessing
global RecursiveStringGeneration_Handle
try
nOut=nargout(RecursiveStringGeneration_Handle);
catch
include_PostProcessing
nOut=nargout(RecursiveStringGeneration_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RecursiveStringGeneration_Handle(varargin{:});
end
