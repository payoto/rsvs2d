function [varargout]=BlueCrystalCheck(varargin)
% include_Utilities
global BlueCrystalCheck_Handle
nOut=nargout(BlueCrystalCheck_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BlueCrystalCheck_Handle(varargin{:});
end
