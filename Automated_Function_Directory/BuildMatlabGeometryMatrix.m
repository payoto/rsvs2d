function [varargout]=BuildMatlabGeometryMatrix(varargin)
global BuildMatlabGeometryMatrix_Handle
nOut=nargout(BuildMatlabGeometryMatrix_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildMatlabGeometryMatrix_Handle(varargin{:});
end
