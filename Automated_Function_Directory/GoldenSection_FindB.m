function [varargout]=GoldenSection_FindB(varargin)
global GoldenSection_FindB_Handle
nOut=nargout(GoldenSection_FindB_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GoldenSection_FindB_Handle(varargin{:});
end
