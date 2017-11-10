function [varargout]=MakeCartesianGridBoundsInactEunit(varargin)
global MakeCartesianGridBoundsInactEunit_Handle
nOut=nargout(MakeCartesianGridBoundsInactEunit_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCartesianGridBoundsInactEunit_Handle(varargin{:});
end
