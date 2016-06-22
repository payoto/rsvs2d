function [varargout]=MakeCartesianGridBounds(varargin)
global MakeCartesianGridBounds_Handle
nOut=nargout(MakeCartesianGridBounds_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCartesianGridBounds_Handle(varargin{:});
end
