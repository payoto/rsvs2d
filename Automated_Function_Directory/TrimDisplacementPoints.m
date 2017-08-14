function [varargout]=TrimDisplacementPoints(varargin)
% include_PostProcessing
global TrimDisplacementPoints_Handle
nOut=nargout(TrimDisplacementPoints_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimDisplacementPoints_Handle(varargin{:});
end
