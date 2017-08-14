function [varargout]=CheckSnakeSensitivityAlgorithm(varargin)
% include_Optimisation
global CheckSnakeSensitivityAlgorithm_Handle
nOut=nargout(CheckSnakeSensitivityAlgorithm_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckSnakeSensitivityAlgorithm_Handle(varargin{:});
end
