function [varargout]=ParallelObjectiveCalc(varargin)
global ParallelObjectiveCalc_Handle
nOut=nargout(ParallelObjectiveCalc_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ParallelObjectiveCalc_Handle(varargin{:});
end
