function [varargout]=CheckResults(varargin)
global CheckResults_Handle
nOut=nargout(CheckResults_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckResults_Handle(varargin{:});
end
