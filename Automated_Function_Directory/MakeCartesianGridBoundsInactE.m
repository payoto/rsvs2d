function [varargout]=MakeCartesianGridBoundsInactE(varargin)
% include_Utilities
global MakeCartesianGridBoundsInactE_Handle
nOut=nargout(MakeCartesianGridBoundsInactE_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCartesianGridBoundsInactE_Handle(varargin{:});
end
