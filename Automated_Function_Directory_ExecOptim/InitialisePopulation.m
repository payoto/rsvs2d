function [varargout]=InitialisePopulation(varargin)
global InitialisePopulation_Handle
nOut=nargout(InitialisePopulation_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialisePopulation_Handle(varargin{:});
end
