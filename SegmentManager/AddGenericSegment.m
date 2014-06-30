function Segment = AddGenericSegment(Segment, newSegmentName, lon1, lat1, lon2, lat2)
%%  AddGenericSegment

Segment.name          = strvcat(Segment.name, newSegmentName);
Segment.lon1          = [Segment.lon1 ; lon1];
Segment.lat1          = [Segment.lat1 ; lat1];
Segment.lon2          = [Segment.lon2 ; lon2];
Segment.lat2          = [Segment.lat2 ; lat2];

Segment.lDep          = [Segment.lDep ; 15];
Segment.lDepSig       = [Segment.lDepSig ; 5];
Segment.lDepTog       = [Segment.lDepTog ; 0];

Segment.dip           = [Segment.dip ; 90];
Segment.dipSig        = [Segment.dipSig ; 1];
Segment.dipTog        = [Segment.dipTog ; 0];

Segment.ssRate        = [Segment.ssRate ; 0];
Segment.ssRateSig     = [Segment.ssRateSig ; 1];
Segment.ssRateTog     = [Segment.ssRateTog ; 0];

Segment.dsRate        = [Segment.dsRate ; 0];
Segment.dsRateSig     = [Segment.dsRateSig ; 1];
Segment.dsRateTog     = [Segment.dsRateTog ; 0];

Segment.tsRate        = [Segment.tsRate ; 0];
Segment.tsRateSig     = [Segment.tsRateSig ; 1];
Segment.tsRateTog     = [Segment.tsRateTog ; 0];

Segment.bDep          = [Segment.bDep ; 15];
Segment.bDepSig       = [Segment.bDepSig ; 5];
Segment.bDepTog       = [Segment.bDepTog ; 0];

Segment.res           = [Segment.res ; 100];
Segment.resOver       = [Segment.resOver ; 0];
Segment.resOther      = [Segment.resOther ; 0];

Segment.patchFile     = [Segment.patchFile ; 0];
Segment.patchTog      = [Segment.patchTog ; 0];
Segment.other3        = [Segment.other3  ; 0];

Segment.patchSlipFile = [Segment.patchSlipFile ; 0];
Segment.patchSlipTog  = [Segment.patchSlipTog  ; 0];
Segment.other6        = [Segment.other6 ; 0];

Segment.other7        = [Segment.other7 ; 0];
Segment.other8        = [Segment.other8 ; 0];
Segment.other9        = [Segment.other9 ; 0];

Segment.other10       = [Segment.other10 ; 0];
Segment.other11       = [Segment.other11 ; 0];
Segment.other12       = [Segment.other12 ; 0];
