function [varargout]=Rosenbrock(varargin)
global Rosenbrock_Handle
nOut=nargout(Rosenbrock_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=Rosenbrock_Handle(varargin{:});
end
