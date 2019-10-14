function [varargout]=ConstantArea_Parabola(varargin)
% include_Optimisation
global ConstantArea_Parabola_Handle
try
nOut=nargout(ConstantArea_Parabola_Handle);
catch
include_Optimisation
nOut=nargout(ConstantArea_Parabola_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConstantArea_Parabola_Handle(varargin{:});
end
