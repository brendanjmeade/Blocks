function Segment = CopySegmentProp(Segment, idx, newSegmentName, lon1, lat1, lon2, lat2)
%%  CopySegmentProp

Segment.name         = strvcat(Segment.name, newSegmentName);
Segment.lon1         = [Segment.lon1 ; lon1];
Segment.lat1         = [Segment.lat1 ; lat1];
Segment.lon2         = [Segment.lon2 ; lon2];
Segment.lat2         = [Segment.lat2 ; lat2];

Segment.lDep         = [Segment.lDep ; Segment.lDep(idx)];
Segment.lDepSig      = [Segment.lDepSig ; Segment.lDepSig(idx)];
Segment.lDepTog      = [Segment.lDepTog ; Segment.lDepTog(idx)];

Segment.dip          = [Segment.dip ; Segment.dip(idx)];
Segment.dipSig       = [Segment.dipSig ; Segment.dipSig(idx)];
Segment.dipTog       = [Segment.dipTog ; Segment.dipTog(idx)];

Segment.ssRate       = [Segment.ssRate ; Segment.ssRate(idx)];
Segment.ssRateSig    = [Segment.ssRateSig ; Segment.ssRateSig(idx)];
Segment.ssRateTog    = [Segment.ssRateTog ; Segment.ssRateTog(idx)];

Segment.dsRate       = [Segment.dsRate ; Segment.dsRate(idx)];
Segment.dsRateSig    = [Segment.dsRateSig ; Segment.dsRateSig(idx)];
Segment.dsRateTog    = [Segment.dsRateTog ; Segment.dsRateTog(idx)];

Segment.tsRate       = [Segment.tsRate ; Segment.tsRate(idx)];
Segment.tsRateSig    = [Segment.tsRateSig ; Segment.tsRateSig(idx)];
Segment.tsRateTog    = [Segment.tsRateTog ; Segment.tsRateTog(idx)];

Segment.bDep         = [Segment.bDep ; Segment.bDep(idx)];
Segment.bDepSig      = [Segment.bDepSig ; Segment.bDepSig(idx)];
Segment.bDepTog      = [Segment.bDepTog ; Segment.bDepTog(idx)];

Segment.res          = [Segment.res ; Segment.res(idx)];
Segment.resOver      = [Segment.resOver ; Segment.resOver(idx)];
Segment.resOther     = [Segment.resOther ; Segment.resOther(idx)];

Segment.other1       = [Segment.other1 ; Segment.other1(idx)];
Segment.other2       = [Segment.other2 ; Segment.other2(idx)];
Segment.other3       = [Segment.other3 ; Segment.other3(idx)];

Segment.other4       = [Segment.other4 ; Segment.other4(idx)];
Segment.other5       = [Segment.other5 ; Segment.other5(idx)];
Segment.other6       = [Segment.other6 ; Segment.other6(idx)];

Segment.other7       = [Segment.other7 ; Segment.other7(idx)];
Segment.other8       = [Segment.other8 ; Segment.other8(idx)];
Segment.other9       = [Segment.other9 ; Segment.other9(idx)];

Segment.other10      = [Segment.other10 ; Segment.other10(idx)];
Segment.other11      = [Segment.other11 ; Segment.other11(idx)];
Segment.other12      = [Segment.other12 ; Segment.other12(idx)];

Segment.patchName    = strvcat(Segment.patchName, Segment.patchName(idx, :));
Segment.patchFlag01  = [Segment.patchFlag01 ; Segment.patchFlag01(idx)];
Segment.patchFlag02  = [Segment.patchFlag02 ; Segment.patchFlag02(idx)];
