function [varargout]=CalculateSnaxelCurvature(varargin)
% include_SnakeParam
global CalculateSnaxelCurvature_Handle
try
nOut=nargout(CalculateSnaxelCurvature_Handle);
catch
include_SnakeParam
nOut=nargout(CalculateSnaxelCurvature_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateSnaxelCurvature_Handle(varargin{:});
end
