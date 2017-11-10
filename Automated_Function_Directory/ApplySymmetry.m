function [varargout]=ApplySymmetry(varargin)
% include_Optimisation
global ApplySymmetry_Handle
nOut=nargout(ApplySymmetry_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ApplySymmetry_Handle(varargin{:});
end
