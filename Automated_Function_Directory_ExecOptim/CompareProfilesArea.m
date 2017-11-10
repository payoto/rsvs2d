function [varargout]=CompareProfilesArea(varargin)
global CompareProfilesArea_Handle
nOut=nargout(CompareProfilesArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesArea_Handle(varargin{:});
end
