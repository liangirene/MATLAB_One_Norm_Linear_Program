function [data, info] = OneNormLP(A,b)

[m, n] = size(A); %Computes problem dimensions

info.run = 'Success';
data.loop = 0;
data.x = NaN * ones(m, 1);
data.obj = NaN;

if (m <= n) %Check dimensions
    info.msg = ('m <= n. Failure due to incorrect dimension of A!')
    info.run = 'Failure'
    return
end

B = 1:n; %Select initial index set with |B| = n

M = A(B,:)^(-1);

opt_flg = 0 %Set optimal flag

while opt_flg == 0
    data.loop = data.loop + 1
    
    if data.loop > nchoosek(m, n)
        info.msg = ('Too many loops')
        info.run = 'Failure'
        return
    end 
    
    if (~isfinite(sum(M(:)))) %Check non-degeneracy assumptions
        info.msg = ("Failure due to degeneracy. Doesn't satisfy non-degeneracy assumption #1")
        info.run = 'Failure'
    end
    
    y = zeros(m, 1);
    B_bar = setdiff(1:m,B);
    x = M*b(B); %Compute x
    h = A*x - b; %Define h. Note that h(B)=0
    h(B_bar) = A(B_bar, :)*x - b(B_bar);
    y(B_bar) = sign(h(B_bar)); %Set y(B_bar) so that y(B_Bar) contains the signs of the h(B_bar) components
    y(B) = -transpose(M) * (transpose(A(B_bar, :)) * y(B_bar)); %Compute y(B)
    
    disp("B:" + B)
    disp("Objective Value:" + ones(1, m) * abs(A*x - b))
    
    if (abs(y(B)) <= 1) %y is feasible solution to dual, so end loop
        opt_flg = 1
        data.obj = ones(1, m) * abs(A*x - b) %optimal objective value
        data.x = x %optimal solution as column vector
        data.B = B
        data.loop = data.loop %no. of iterations to solve for optimal solution
        return
    end
        
    y_index = find(abs(y(B)) > 1, 1); %Compute index where abs(y_s) > 1
    s = y_index; %Choose s in B such that abs(y(B))>1. j_s leaves B
    
    t = zeros(m, 1);
    t(B_bar) = -(sign(y(B(s))))*(y(B_bar)).*(A(B_bar,:)* M(:, s)); %Compute t(B_Bar)

    mask = find(t>0); %filter so that t_j > 0  
    [min_val_r, r] = min(abs(h(mask))./(t(mask))); %Compute r. r joins B.
    r = mask(r);
    
    if h(mask) == 0 %Check that there does not exist an index set B with more than n indexes such that A(B,:)*x = b(B)
        info.msg = ("Failure due to degeneracy. Doesn't satisfy non-degeneracy assumption #2")
        info.run = 'Failure'
    end 
    
    B(s) = r; %Update index set. New active basis, kicking out j_s, adding r
    
    theta = transpose(A(r,:) * M); %compute theta^T
    M(:, s) = (1/theta(s)) * M(:,s); %Update inverse matrix
    
    [c, d] = size(B);
    B_indices = 1:d;
    j_not_s = B_indices(B_indices~=s); %for j =/= s
    for j = j_not_s
        M(:,j) = M(:,j) - (theta(j) * M(:,s));
    end
    
end
end
