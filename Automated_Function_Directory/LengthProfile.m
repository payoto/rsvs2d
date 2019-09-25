function [varargout]=LengthProfile(varargin)
% include_Utilities
global LengthProfile_Handle
try
nOut=nargout(LengthProfile_Handle);
catch
include_Utilities
nOut=nargout(LengthProfile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LengthProfile_Handle(varargin{:});
end
