function [varargout]=PointGeneration(varargin)
% include_SnakeParam
global PointGeneration_Handle
try
nOut=nargout(PointGeneration_Handle);
catch
include_SnakeParam
nOut=nargout(PointGeneration_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PointGeneration_Handle(varargin{:});
end
