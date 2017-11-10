function [varargout]=CompareProfilesDistance2(varargin)
global CompareProfilesDistance2_Handle
nOut=nargout(CompareProfilesDistance2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareProfilesDistance2_Handle(varargin{:});
end
