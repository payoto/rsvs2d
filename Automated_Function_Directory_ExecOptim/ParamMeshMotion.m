function [varargout]=ParamMeshMotion(varargin)
global ParamMeshMotion_Handle
nOut=nargout(ParamMeshMotion_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ParamMeshMotion_Handle(varargin{:});
end
