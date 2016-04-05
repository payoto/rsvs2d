function [varargout]=ConcatenateAutomaticComments(varargin)
global ConcatenateAutomaticComments_Handle
nOut=nargout(ConcatenateAutomaticComments_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConcatenateAutomaticComments_Handle(varargin{:});
end
