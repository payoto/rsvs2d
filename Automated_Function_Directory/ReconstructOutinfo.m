function [varargout]=ReconstructOutinfo(varargin)
% include_PostProcessing
global ReconstructOutinfo_Handle
try
nOut=nargout(ReconstructOutinfo_Handle);
catch
include_PostProcessing
nOut=nargout(ReconstructOutinfo_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReconstructOutinfo_Handle(varargin{:});
end
