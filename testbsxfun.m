v = [1,  2,  3];

A = [8,  1,  6
     3,  5,  7
     4,  9,  2];
A=repmat(A,1000000,1);
 
ii=1;

v = [1,  2,  3];
t(ii,1)=now;
B = bsxfun(@plus, A, v);
t(ii,2)=now;
ii=ii+1;

t(ii,1)=now;
% B = zeros(size(A));
% for row = 1:size(A,1)
%     B(row,:) = A(row,:) + v;
% end
t(ii,2)=now;
ii=ii+1;

v = [1,  2,  3];


t(ii,1)=now;
 
B = A + repmat(v,size(A,1),1); 
t(ii,2)=now;
ii=ii+1;




Dt=t(:,2)-t(:,1);