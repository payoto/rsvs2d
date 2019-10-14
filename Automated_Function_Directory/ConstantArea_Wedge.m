function [varargout]=ConstantArea_Wedge(varargin)
% include_Optimisation
global ConstantArea_Wedge_Handle
try
nOut=nargout(ConstantArea_Wedge_Handle);
catch
include_Optimisation
nOut=nargout(ConstantArea_Wedge_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConstantArea_Wedge_Handle(varargin{:});
end
