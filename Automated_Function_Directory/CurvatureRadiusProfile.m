function [varargout]=CurvatureRadiusProfile(varargin)
% include_Utilities
global CurvatureRadiusProfile_Handle
try
nOut=nargout(CurvatureRadiusProfile_Handle);
catch
include_Utilities
nOut=nargout(CurvatureRadiusProfile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CurvatureRadiusProfile_Handle(varargin{:});
end
