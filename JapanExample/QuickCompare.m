function QuickCompare(dir1, dir2)
% Quick comparison of Blocks output files for testing purposes.

% Compare segment files
Seg1 = ReadSegmentStruct(strcat(dir1, '/Mod.segment'));
Seg2 = ReadSegmentStruct(strcat(dir2, '/Mod.segment'));

% Are there the same number of segments
if numel(Seg1.lon1) == numel(Seg2.lon1)
    fprintf(1, 'Number segments match.  Have a nice day.\n')
else
    fprintf(1, 'Number segments do not match.\n')
    fprintf(1, '%s has %d segments.\n', dir1, numel(Seg1.lon1));
    fprintf(1, '%s has %d segments.\n', dir2, numel(Seg2.lon1));
end


for i = 1:numel(Seg1.lon1)
    if Seg1.ssRate(i) ~= Seg2.ssRate(i)
       fprintf(1, '%s has different ssRate.\n', deblank(Seg1.name(i,:)));
       fprintf(1, '%s ssRate %6.3f.\n', dir1, Seg1.ssRate(i));
       fprintf(1, '%s ssRate %6.3f.\n', dir2, Seg2.ssRate(i));
    end

    if Seg1.dsRate(i) ~= Seg2.dsRate(i)
       fprintf(1, '%s has different dsRate.\n', deblank(Seg1.name(i,:)));
       fprintf(1, '%s dsRate %6.3f.\n', dir1, Seg1.dsRate(i));
       fprintf(1, '%s dsRate %6.3f.\n', dir2, Seg2.dsRate(i));
    end

    if Seg1.tsRate(i) ~= Seg2.tsRate(i)
       fprintf(1, '%s has different ssRate.\n', deblank(Seg1.name(i,:)));
       fprintf(1, '%s tsRate %6.3f.\n', dir1, Seg1.tsRate(i));
       fprintf(1, '%s tsRate %6.3f.\n', dir2, Seg2.tsRate(i));
    end
end
