function [varargout]=FindObjNum(varargin)
global FindObjNum_Handle
nOut=nargout(FindObjNum_Handle);
[varargout{1:nOut}]=FindObjNum_Handle(varargin{:});
end
