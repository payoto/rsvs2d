function [varargout]=DisplacementsOutput(varargin)
% include_PostProcessing
global DisplacementsOutput_Handle
try
nOut=nargout(DisplacementsOutput_Handle);
catch
include_PostProcessing
nOut=nargout(DisplacementsOutput_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DisplacementsOutput_Handle(varargin{:});
end
