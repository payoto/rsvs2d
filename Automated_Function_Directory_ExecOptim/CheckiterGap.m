function [varargout]=CheckiterGap(varargin)
global CheckiterGap_Handle
nOut=nargout(CheckiterGap_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckiterGap_Handle(varargin{:});
end
