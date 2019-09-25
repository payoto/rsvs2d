function [varargout]=HandleNonFillVar(varargin)
% include_Optimisation
global HandleNonFillVar_Handle
try
nOut=nargout(HandleNonFillVar_Handle);
catch
include_Optimisation
nOut=nargout(HandleNonFillVar_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=HandleNonFillVar_Handle(varargin{:});
end
