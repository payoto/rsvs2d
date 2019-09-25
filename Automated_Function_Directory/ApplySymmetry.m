function [varargout]=ApplySymmetry(varargin)
% include_Optimisation
global ApplySymmetry_Handle
try
nOut=nargout(ApplySymmetry_Handle);
catch
include_Optimisation
nOut=nargout(ApplySymmetry_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ApplySymmetry_Handle(varargin{:});
end
