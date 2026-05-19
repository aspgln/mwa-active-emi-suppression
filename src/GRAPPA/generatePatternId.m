function [patternId, patternList, maxPatternId] = generatePatternId(kData, kSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the number of undersampling patterns in the
% grappa data.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sx,sy,nCoil] = size(kData);
kData_tmp = zeroPadding(kData,[sx+kSize(1)-1, sy+kSize(2)-1,nCoil]);

pat_num_limit = 50;
patternId = zeros(sx,sy);
patternList = zeros(kSize(1),kSize(2),nCoil,pat_num_limit);
maxPatternId = 1;
for y = 1:sy
    for x = 1:sx
        tmp = kData_tmp(x:x+kSize(1)-1,y:y+kSize(2)-1,:);
        pat = abs(tmp)>0;
        count_num = 1;
        find_pat = 0;
        while find_pat == 0
            if sum(patternList(:,:,:,count_num)) == 0
                patternList(:,:,:,count_num) = pat;
                find_pat = 1;
                patternId(x,y) = count_num;
                if count_num > maxPatternId
                    maxPatternId = count_num;
                end
            else
                if isequal(pat,patternList(:,:,:,count_num))
                    patternId(x,y) = count_num;
                    if count_num > maxPatternId
                        maxPatternId = count_num;
                    end
                    find_pat = 1;
                else
                    count_num = count_num + 1;
                end

            end
        end
    end
end