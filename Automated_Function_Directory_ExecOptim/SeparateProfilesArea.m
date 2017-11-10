function [varargout]=SeparateProfilesArea(varargin)
global SeparateProfilesArea_Handle
nOut=nargout(SeparateProfilesArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SeparateProfilesArea_Handle(varargin{:});
end
