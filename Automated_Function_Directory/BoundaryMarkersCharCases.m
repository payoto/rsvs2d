function [varargout]=BoundaryMarkersCharCases(varargin)
% include_PostProcessing
global BoundaryMarkersCharCases_Handle
try
nOut=nargout(BoundaryMarkersCharCases_Handle);
catch
include_PostProcessing
nOut=nargout(BoundaryMarkersCharCases_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BoundaryMarkersCharCases_Handle(varargin{:});
end
