function [varargout]=RemoveIdenticalEntries(varargin)
global RemoveIdenticalEntries_Handle
nOut=nargout(RemoveIdenticalEntries_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalEntries_Handle(varargin{:});
end
