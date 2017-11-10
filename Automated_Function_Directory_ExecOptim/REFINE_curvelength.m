function [varargout]=REFINE_curvelength(varargin)
global REFINE_curvelength_Handle
nOut=nargout(REFINE_curvelength_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_curvelength_Handle(varargin{:});
end
