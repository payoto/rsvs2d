function [varargout]=BoundaryInput(varargin)
global BoundaryInput_Handle
nOut=nargout(BoundaryInput_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BoundaryInput_Handle(varargin{:});
end
