function [varargout]=InitialiseObjective(varargin)
global InitialiseObjective_Handle
nOut=nargout(InitialiseObjective_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseObjective_Handle(varargin{:});
end
