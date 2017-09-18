function [varargout]=MakeBoundsOuterLayer(varargin)
% include_Utilities
global MakeBoundsOuterLayer_Handle
nOut=nargout(MakeBoundsOuterLayer_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeBoundsOuterLayer_Handle(varargin{:});
end
