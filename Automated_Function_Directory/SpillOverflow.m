function [varargout]=SpillOverflow(varargin)
global SpillOverflow_Handle
nOut=nargout(SpillOverflow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SpillOverflow_Handle(varargin{:});
end
