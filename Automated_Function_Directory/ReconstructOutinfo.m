function [varargout]=ReconstructOutinfo(varargin)
% include_PostProcessing
global ReconstructOutinfo_Handle
nOut=nargout(ReconstructOutinfo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReconstructOutinfo_Handle(varargin{:});
end
