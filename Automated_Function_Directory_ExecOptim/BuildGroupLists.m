function [varargout]=BuildGroupLists(varargin)
global BuildGroupLists_Handle
nOut=nargout(BuildGroupLists_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildGroupLists_Handle(varargin{:});
end
