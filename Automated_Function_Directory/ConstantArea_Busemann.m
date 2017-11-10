function [varargout]=ConstantArea_Busemann(varargin)
% include_Optimisation
global ConstantArea_Busemann_Handle
nOut=nargout(ConstantArea_Busemann_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConstantArea_Busemann_Handle(varargin{:});
end
