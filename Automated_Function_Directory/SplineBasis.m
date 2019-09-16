function [varargout]=SplineBasis(varargin)
% include_NURBSEngine
global SplineBasis_Handle
nOut=nargout(SplineBasis_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SplineBasis_Handle(varargin{:});
end
