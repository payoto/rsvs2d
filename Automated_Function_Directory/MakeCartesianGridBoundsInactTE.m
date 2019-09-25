function [varargout]=MakeCartesianGridBoundsInactTE(varargin)
% include_Utilities
global MakeCartesianGridBoundsInactTE_Handle
try
nOut=nargout(MakeCartesianGridBoundsInactTE_Handle);
catch
include_Utilities
nOut=nargout(MakeCartesianGridBoundsInactTE_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCartesianGridBoundsInactTE_Handle(varargin{:});
end
