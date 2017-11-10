function [varargout]=BoundaryMarkersCharCases(varargin)
% include_PostProcessing
global BoundaryMarkersCharCases_Handle
nOut=nargout(BoundaryMarkersCharCases_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BoundaryMarkersCharCases_Handle(varargin{:});
end
