function [varargout]=MakeCartesianGridBoundsInactTE(varargin)
global MakeCartesianGridBoundsInactTE_Handle
nOut=nargout(MakeCartesianGridBoundsInactTE_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCartesianGridBoundsInactTE_Handle(varargin{:});
end
