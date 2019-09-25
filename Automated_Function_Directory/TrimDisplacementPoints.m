function [varargout]=TrimDisplacementPoints(varargin)
% include_PostProcessing
global TrimDisplacementPoints_Handle
try
nOut=nargout(TrimDisplacementPoints_Handle);
catch
include_PostProcessing
nOut=nargout(TrimDisplacementPoints_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimDisplacementPoints_Handle(varargin{:});
end
