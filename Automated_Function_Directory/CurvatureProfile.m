function [varargout]=CurvatureProfile(varargin)
% include_Utilities
global CurvatureProfile_Handle
try
nOut=nargout(CurvatureProfile_Handle);
catch
include_Utilities
nOut=nargout(CurvatureProfile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CurvatureProfile_Handle(varargin{:});
end
