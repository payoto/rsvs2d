function [varargout]=MakeBoundsOuterLayer(varargin)
% include_Utilities
global MakeBoundsOuterLayer_Handle
try
nOut=nargout(MakeBoundsOuterLayer_Handle);
catch
include_Utilities
nOut=nargout(MakeBoundsOuterLayer_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeBoundsOuterLayer_Handle(varargin{:});
end
