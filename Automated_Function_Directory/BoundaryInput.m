function [varargout]=BoundaryInput(varargin)
% include_PostProcessing
global BoundaryInput_Handle
try
nOut=nargout(BoundaryInput_Handle);
catch
include_PostProcessing
nOut=nargout(BoundaryInput_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BoundaryInput_Handle(varargin{:});
end
