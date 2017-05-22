function [varargout]=CompareProfilesDistance(varargin)
global CompareProfilesDistance_Handle
nOut=nargout(CompareProfilesDistance_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesDistance_Handle(varargin{:});
end
