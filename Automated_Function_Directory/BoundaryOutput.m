function [varargout]=BoundaryOutput(varargin)
% include_PostProcessing
global BoundaryOutput_Handle
nOut=nargout(BoundaryOutput_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BoundaryOutput_Handle(varargin{:});
end
