function [varargout]=DisplacementsOutput(varargin)
% include_PostProcessing
global DisplacementsOutput_Handle
nOut=nargout(DisplacementsOutput_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DisplacementsOutput_Handle(varargin{:});
end
