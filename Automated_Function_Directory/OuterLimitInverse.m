function [varargout]=OuterLimitInverse(varargin)
% include_Optimisation
global OuterLimitInverse_Handle
nOut=nargout(OuterLimitInverse_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OuterLimitInverse_Handle(varargin{:});
end
