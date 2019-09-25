function [varargout]=savefig(varargin)
% include_PostProcessing
global savefig_Handle
try
nOut=nargout(savefig_Handle);
catch
include_PostProcessing
nOut=nargout(savefig_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=savefig_Handle(varargin{:});
end
