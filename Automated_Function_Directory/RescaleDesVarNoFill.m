function [varargout]=RescaleDesVarNoFill(varargin)
global RescaleDesVarNoFill_Handle
nOut=nargout(RescaleDesVarNoFill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RescaleDesVarNoFill_Handle(varargin{:});
end
