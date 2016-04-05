function [varargout]=LengthProfile(varargin)
global LengthProfile_Handle
nOut=nargout(LengthProfile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LengthProfile_Handle(varargin{:});
end
