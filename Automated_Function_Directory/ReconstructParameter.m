function [varargout]=ReconstructParameter(varargin)
% include_PostProcessing
global ReconstructParameter_Handle
try
nOut=nargout(ReconstructParameter_Handle);
catch
include_PostProcessing
nOut=nargout(ReconstructParameter_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReconstructParameter_Handle(varargin{:});
end
