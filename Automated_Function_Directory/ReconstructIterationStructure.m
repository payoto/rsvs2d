function [varargout]=ReconstructIterationStructure(varargin)
% include_PostProcessing
global ReconstructIterationStructure_Handle
try
nOut=nargout(ReconstructIterationStructure_Handle);
catch
include_PostProcessing
nOut=nargout(ReconstructIterationStructure_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReconstructIterationStructure_Handle(varargin{:});
end
