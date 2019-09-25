function [varargout]=BlueCrystalCheck(varargin)
% include_Utilities
global BlueCrystalCheck_Handle
try
nOut=nargout(BlueCrystalCheck_Handle);
catch
include_Utilities
nOut=nargout(BlueCrystalCheck_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BlueCrystalCheck_Handle(varargin{:});
end
