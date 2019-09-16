function [varargout]=SerialisedParallelObjectiveCalc(varargin)
% ExecuteOptimisation
global SerialisedParallelObjectiveCalc_Handle
nOut=nargout(SerialisedParallelObjectiveCalc_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SerialisedParallelObjectiveCalc_Handle(varargin{:});
end
