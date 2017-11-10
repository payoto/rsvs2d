function [varargout]=ReconstructParameter(varargin)
% include_PostProcessing
global ReconstructParameter_Handle
nOut=nargout(ReconstructParameter_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReconstructParameter_Handle(varargin{:});
end
