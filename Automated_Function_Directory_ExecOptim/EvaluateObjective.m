function [varargout]=EvaluateObjective(varargin)
global EvaluateObjective_Handle
nOut=nargout(EvaluateObjective_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EvaluateObjective_Handle(varargin{:});
end
