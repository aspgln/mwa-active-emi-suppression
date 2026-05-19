function [order,order_rev] = permute2RefDim(dimStr,refDimStr)

Nstr = length(dimStr);
order = zeros(1,Nstr);

count = Nstr+1;
for i = 1:length(refDimStr)
    strLoc = strfind(dimStr,refDimStr(i));
    if ~isempty(strLoc)
        order(i) = strLoc;
    else
        order(i) = count;
        count = count + 1;
    end
end


[~, order_rev] = sort(order,'ascend');
