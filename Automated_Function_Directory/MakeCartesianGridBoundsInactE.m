function [varargout]=MakeCartesianGridBoundsInactE(varargin)
% include_Utilities
global MakeCartesianGridBoundsInactE_Handle
try
nOut=nargout(MakeCartesianGridBoundsInactE_Handle);
catch
include_Utilities
nOut=nargout(MakeCartesianGridBoundsInactE_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCartesianGridBoundsInactE_Handle(varargin{:});
end
