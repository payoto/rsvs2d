function [varargout]=CurvatureRadiusProfile(varargin)
% include_Utilities
global CurvatureRadiusProfile_Handle
nOut=nargout(CurvatureRadiusProfile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CurvatureRadiusProfile_Handle(varargin{:});
end
