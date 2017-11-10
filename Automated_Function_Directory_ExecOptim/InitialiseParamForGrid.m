function [varargout]=InitialiseParamForGrid(varargin)
global InitialiseParamForGrid_Handle
nOut=nargout(InitialiseParamForGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseParamForGrid_Handle(varargin{:});
end
