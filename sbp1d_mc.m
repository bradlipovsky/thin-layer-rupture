function mats = sbp1d_mc(N,h,order)
% creates compatible 1D SBP factors for creation of 2D operators


if N == 1,
  if order == 2,
    mats.H = speye(N,N);
    mats.Hinv = sparse(N,N);
    mats.BS = sparse(N,N);
    mats.D1int = sparse(N,N);
    mats.D1 = sparse(N,N);
    mats.D2 = sparse(N,N);
  elseif order == 4
    mats.H = speye(N,N);
    mats.Hinv = sparse(N,N);
    mats.BS = sparse(N,N);
    mats.D1int = sparse(N,N);
    mats.D1 = sparse(N,N);
    mats.D3 = sparse(N,N);
    mats.D4 = sparse(N,N);
    mats.C3 = sparse(N,N);
    mats.C4 = sparse(N,N);
  end
  return
end

switch order
  case(2)
    mats.H = sparse(1:N,1:N,[0.5,ones(1,N-2),0.5]).*h;
    mats.Hinv = sparse(1:N,1:N,[2,ones(1,N-2),2])./h;
    
    % 1st row is -1*p666 of Mattsson 2010
    mats.BS = sparse([1 1 1],[1 2 3],[1.5 -2 0.5],N,N)./h +...
      sparse([N N N],[N-2 N-1 N],[0.5 -2 1.5],N,N)./h;
    
    % interior 1st derivative
    mats.D1int = sparse(2:N-1,1:N-2, -0.5*ones(1,N-2),N,N) + ...
      sparse(2:N-1,3:N, 0.5*ones(1,N-2),N,N);
    mats.D1int(1,[1 2]) = [-1 1];
    mats.D1int(end,[end-1 end]) = [-1 1];
    mats.D1int = mats.D1int./h;
    
    % 1st derivative with transition to one-sided at boundaries
    mats.D1 = mats.D1int;
    mats.D1(1,:) = 0.*mats.D1(1,:); mats.D1(1,1:3) = -mats.BS(1,1:3);
    mats.D1(end,:) = 0.*mats.D1(end,:); mats.D1(end,end-3:end) = mats.BS(end,end-3:end);
    
    % 2nd derivative
    mats.D2 = sparse(2:N-1,2:N-1,-2.*ones(1,N-2),N,N) + ...
      sparse(2:N-1,3:N,ones(1,N-2),N,N) + ...
      sparse(2:N-1,1:N-2,ones(1,N-2),N,N);
    mats.D2(1,1:3) = [1 -2 1];
    mats.D2(N,N-2:N) = [1 -2 1];
    mats.D2 = mats.D2;
    
  case(4)
    vec = [17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0];
    mats.H = sparse(1:N,1:N,[vec,ones(1,N-8),fliplr(vec)]).*h;
    mats.Hinv = inv(mats.H);
    
    
    % row 1: -1* p666 of Mattsson 2010
    mats.BS = sparse(N,N);
    vec = [11.0/6.0, -3.0, 1.5, -1/3];
    mats.BS(1,1:4) = vec;
    mats.BS(N,N-3:N) = fliplr(vec);
    mats.BS = mats.BS./h;
    
    % interior stencil for 1st derivative
    rowV = 5:N-4;
    colV = 3:N-6;
    mats.D1int = sparse(rowV,colV,1/12*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+1,-2/3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+3,2/3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+4,-1/12*ones(size(rowV)),N,N);
    
    mats.D1int(1:4,1:6) = [-24/17,59/34,-4/17,-3/34,0,0; -1/2,0,1/2,0,0,0;...
      4/43,-59/86,0,59/86,-4/43,0; 3/98,0,-59/98,0,32/49,-4/49];
    mats.D1int(N-3:N,N-5:N) = rot90( -mats.D1int(1:4,1:6),2);
    mats.D1int = mats.D1int./h;
    
    mats.D1= mats.D1int;
    mats.D1(1,:) = 0.*mats.D1(1,:); mats.D1(1,:) = -mats.BS(1,:);
    mats.D1(end,:) = 0.*mats.D1(end,:); mats.D1(N,:) = mats.BS(N,:);
    
    
    % repeating interior stencil [-1, 3, -3, 1]
    rowV = [1,2,4:N-3,N-1,N];
    colV = [1,1,3:N-4,N-3,N-3];
    mats.D3 = sparse(rowV,colV,-1*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+1,3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+2,-3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+3,1*ones(size(rowV)),N,N);
    % rows 3 and N-3
    vec = [-185893.0/301051.0, 79000249461.0/54642863857.0, -33235054191.0/54642863857.0,...
      -36887526683.0/54642863857.0, 26183621850.0/54642863857.0, -4386.0/181507.0];
    mats.D3(3,1:6) = vec;
    mats.D3(N-2,N-5:N) = -fliplr(vec);
    
    rowV = 1:N;
    colV = [1,1,1,2:N-5,N-4,N-4,N-4];
    mats.D4 = sparse(rowV,colV,1*ones(1,N),N,N)...
      + sparse(rowV,colV+1,-4*ones(1,N),N,N)...
      + sparse(rowV,colV+2,6*ones(1,N),N,N)...
      + sparse(rowV,colV+3,-4*ones(1,N),N,N)...
      + sparse(rowV,colV+4,1*ones(1,N),N,N);
    
    vec1 = [0, 0, 163928591571.0/53268010936.0, 189284.0/185893.0];
    vec2 = [189284.0/185893.0, 0, 163928591571.0/53268010936.0, 0, 0];
    mats.C3 = sparse(1:N,1:N,[vec1,ones(1,N-9),vec2],N,N);
    
    vec = [0, 0, 1644330.0/301051.0, 156114.0/181507.0];
    mats.C4 = sparse(1:N,1:N,[vec,ones(1,N-8),fliplr(vec)],N,N);
    
  otherwise
    display('Error: order can only be 2 or 4.')
    return
end