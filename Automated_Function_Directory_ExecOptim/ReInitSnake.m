function [varargout]=ReInitSnake(varargin)
% ExecuteOptimisation
global ReInitSnake_Handle
nOut=nargout(ReInitSnake_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReInitSnake_Handle(varargin{:});
end
