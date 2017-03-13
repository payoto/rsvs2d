function [varargout]=REFINE_contlengthnorm(varargin)
global REFINE_contlengthnorm_Handle
nOut=nargout(REFINE_contlengthnorm_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contlengthnorm_Handle(varargin{:});
end
