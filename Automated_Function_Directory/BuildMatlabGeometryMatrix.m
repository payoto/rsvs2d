function [varargout]=BuildMatlabGeometryMatrix(varargin)
% include_PostProcessing
global BuildMatlabGeometryMatrix_Handle
try
nOut=nargout(BuildMatlabGeometryMatrix_Handle);
catch
include_PostProcessing
nOut=nargout(BuildMatlabGeometryMatrix_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildMatlabGeometryMatrix_Handle(varargin{:});
end
