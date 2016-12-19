function [varargout]=BuildCellRatios(varargin)
global BuildCellRatios_Handle
nOut=nargout(BuildCellRatios_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildCellRatios_Handle(varargin{:});
end
