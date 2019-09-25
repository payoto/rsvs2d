function [varargout]=ConcatenateAutomaticComments(varargin)
% include_PostProcessing
global ConcatenateAutomaticComments_Handle
try
nOut=nargout(ConcatenateAutomaticComments_Handle);
catch
include_PostProcessing
nOut=nargout(ConcatenateAutomaticComments_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConcatenateAutomaticComments_Handle(varargin{:});
end
