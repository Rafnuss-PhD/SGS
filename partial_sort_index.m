% ------------------------------------------------------------------------
% returns the index of the largest N values in the array
function max_N_index = partial_sort_index(x,N)

if(length(x)<N); N=length(x); end
for ii=1:N
    [~, max_N_index(ii)] = max(x); %#ok<AGROW>
    x(max_N_index(ii)) = -inf;
end