function [varargout]=LocalConstraintExtraction(varargin)
global LocalConstraintExtraction_Handle
nOut=nargout(LocalConstraintExtraction_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LocalConstraintExtraction_Handle(varargin{:});
end
