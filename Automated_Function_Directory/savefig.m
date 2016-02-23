function [varargout]=savefig(varargin)
global savefig_Handle
nOut=nargout(savefig_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=savefig_Handle(varargin{:});
end
