function [varargout]=BoundaryOutput(varargin)
% include_PostProcessing
global BoundaryOutput_Handle
try
nOut=nargout(BoundaryOutput_Handle);
catch
include_PostProcessing
nOut=nargout(BoundaryOutput_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BoundaryOutput_Handle(varargin{:});
end
