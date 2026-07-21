function ANS = repop(varargin)
% REPOP computes repeated operations
%    ANS = repop(varargin)
%
%    E.g. repop(x,y,'+',z,'*') returns (x+y)*z

Narg = nargin;
if Narg<3, error('Arguments should be three or more.'), end

Stack = cell(Narg,1); j = 0;

for k = 1:Narg
  arg = varargin{k};
  if isnumeric(arg)
    j = j+1;
    Stack{j} = arg;
  elseif ischar(arg)
    j = j-1;
    Stack{j} = op(Stack{j},Stack{j+1},arg);
  end
end

if j ~= 1, error('Stack Error'), end
ANS = Stack{1};

%%%% Sub Function for REPOP
function ANS = op(A,B,op)

a = size(A);
b = size(B);

dim_a = length(a);
dim_b = length(b);
dim_max = max(dim_a,dim_b);
if dim_max > 2
  a = [a ones(1,dim_max-dim_a)];
  b = [b ones(1,dim_max-dim_b)];
end

if max(a)~=1 | max(b)~=1
  Min = min([a; b]);
  Max = max([a; b]);

  if sum(Min>1 & Min~=Max) > 0, error('Matrix Size Mismatch!'), end
  A = repmat(A,Max./a);
  B = repmat(B,Max./b);
end

switch op
 case '+'
  ANS = A + B;
 case '-'
  ANS = A - B;
 case '*'
  ANS = A .* B;
 case '/'
  ANS = A ./ B;
 case '^'
  ANS = A .^ B;
 case '<'
  ANS = (A < B);
 case '>'
  ANS = (A > B);
 case '<='
  ANS = (A <= B);
 case '>='
  ANS = (A >= B);
 otherwise
  error('Undefined Operator!');
end
