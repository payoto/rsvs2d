function [varargout]=SpillOverflowVarHandling(varargin)
global SpillOverflowVarHandling_Handle
nOut=nargout(SpillOverflowVarHandling_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SpillOverflowVarHandling_Handle(varargin{:});
end
