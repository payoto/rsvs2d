function [varargout]=REFINE_contlength(varargin)
global REFINE_contlength_Handle
nOut=nargout(REFINE_contlength_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contlength_Handle(varargin{:});
end
