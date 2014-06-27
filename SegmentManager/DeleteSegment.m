function Segment = DeleteSegment(Segment, idx)
%%  Delete a segment or a bunch

Segment.name(idx, :)                                                      = [];
Segment.lon1(idx)                                                         = [];
Segment.lat1(idx)                                                         = [];
Segment.lon2(idx)                                                         = [];
Segment.lat2(idx)                                                         = [];

Segment.lDep(idx)                                                         = [];
Segment.lDepSig(idx)                                                      = [];
Segment.lDepTog(idx)                                                      = [];

Segment.bDep(idx)                                                         = [];
Segment.bDepSig(idx)                                                      = [];
Segment.bDepTog(idx)                                                      = [];

Segment.dip(idx)                                                          = [];
Segment.dipSig(idx)                                                       = [];
Segment.dipTog(idx)                                                       = [];

Segment.ssRate(idx)                                                       = [];
Segment.ssRateSig(idx)                                                    = [];
Segment.ssRateTog(idx)                                                    = [];

Segment.dsRate(idx)                                                       = [];
Segment.dsRateSig(idx)                                                    = [];
Segment.dsRateTog(idx)                                                    = [];

Segment.tsRate(idx)                                                       = [];
Segment.tsRateSig(idx)                                                    = [];
Segment.tsRateTog(idx)                                                    = [];

Segment.res(idx)                                                          = [];
Segment.resOver(idx)                                                      = [];
Segment.resOther(idx)                                                     = [];

Segment.other1(idx)                                                       = [];
Segment.other2(idx)                                                       = [];
Segment.other3(idx)                                                       = [];
Segment.other4(idx)                                                       = [];
Segment.other5(idx)                                                       = [];
Segment.other6(idx)                                                       = [];
Segment.other7(idx)                                                       = [];
Segment.other8(idx)                                                       = [];
Segment.other9(idx)                                                       = [];
Segment.other10(idx)                                                      = [];
Segment.other11(idx)                                                      = [];
Segment.other12(idx)                                                      = [];

Segment.patchName(idx, :)                                                 = [];
Segment.patchFlag01(idx)                                                  = [];
Segment.patchFlag02(idx)                                                  = [];