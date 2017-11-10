function [varargout]=CurvatureProfile(varargin)
% include_Utilities
global CurvatureProfile_Handle
nOut=nargout(CurvatureProfile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CurvatureProfile_Handle(varargin{:});
end
