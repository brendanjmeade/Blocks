function WriteSegment(segfilename, Segment)
%%  WriteSegment.m

%%  Open file stream
filestream                          = fopen(segfilename, 'w');

%%  Write header
fprintf(filestream, 'Name\n');
fprintf(filestream, 'st_long  st_lat   end_long  end_lat\n');
fprintf(filestream, 'lock_dep ld_sig   ld_toggle\n');
fprintf(filestream, 'dip      dip_sig  dip_toggle\n');
fprintf(filestream, 'ss_rate  ss_sig   ss_tog\n');
fprintf(filestream, 'ds_rate  ds_sig   ds_tog\n');
fprintf(filestream, 'ts_rate  ts_sig   ts_tog\n');
fprintf(filestream, 'bur_dep  bd_sig   bd_tog\n');
fprintf(filestream, 'fres     fres_ov  fres_other\n');
fprintf(filestream, 'other1   other2   other3\n');
fprintf(filestream, 'other4   other5   other6\n');
fprintf(filestream, 'other7   other8   other9\n');
fprintf(filestream, 'other10  other11  other12\n');

for cnt = 1 : length(Segment.lon1)
   fprintf(filestream, '%s\n', Segment.name(cnt, :));
   fprintf(filestream, '%3.3f   %3.3f   %3.3f  %3.3f\n',   Segment.lon1(cnt),     Segment.lat1(cnt),        Segment.lon2(cnt),        Segment.lat2(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.lDep(cnt),     Segment.lDepSig(cnt),     Segment.lDepTog(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.dip(cnt),      Segment.dipSig(cnt),      Segment.dipTog(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.ssRate(cnt),   Segment.ssRateSig(cnt),   Segment.ssRateTog(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.dsRate(cnt),   Segment.dsRateSig(cnt),   Segment.dsRateTog(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.tsRate(cnt),   Segment.tsRateSig(cnt),   Segment.tsRateTog(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.bDep(cnt),     Segment.bDepSig(cnt),     Segment.bDepTog(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.res(cnt),      Segment.resOver(cnt),     Segment.resOther(cnt));   
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.other1(cnt),   Segment.other2(cnt),      Segment.other3(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.other4(cnt),   Segment.other5(cnt),      Segment.other6(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.other7(cnt),   Segment.other8(cnt),      Segment.other9(cnt));
   fprintf(filestream, '%3.1f   %3.1f  %3.1f\n',           Segment.other10(cnt),  Segment.other11(cnt),     Segment.other12(cnt));
end

%%  Close file
fclose(filestream);
