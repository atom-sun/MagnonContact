On[Assert];

Clear["Global`*"];
Remove["Global`*"];

path0 = DirectoryName[$InputFileName];
path1 = NotebookDirectory[];
If[! DirectoryQ[path0], SetDirectory[path1], SetDirectory[path0]];
Print[Directory[]]


Plot[Sin[x], {x, -10, 10}]


(* Setup plot basics. *)
(*colors=System`PlotThemeDump`$ThemeDefaultGradient*)
Table[{i, ColorData[i] /@ Range[10]}, {i, 115}];
colors0 = ColorData[110] /@ Range[10];
markers = {circle, uptriangle, diamond, square, downtriangle} =
   Charting`CommonDump`GraphicsOpenPlotMarkers[][[;; 5]];
markers = Join[markers, {\[FivePointedStar], \[SixPointedStar]}];
aspratio = 4/5;
imgsize = 360;
frame1[legend_] := Framed[legend,
   Background -> White,
   FrameStyle -> Black,
   FrameMargins -> 0,
   RoundingRadius -> 10
   ];


(* Part I. Demonstrate correlation functions. *)
(* load MPS data. *)
(* model parameters. *)
nSites = 200;
alpha = 1.9;
Jz = 1.9;
z := alpha - 1;

numMagnons = 3;

f = FileNameJoin[{Directory[], "data", 
   StringTemplate["magnon-num_magnons``-n``-a``-Jz``-bc-MPS.h5"][
    numMagnons, nSites, alpha, Jz]}]

zz = Import[f, "zz"];
pm = Import[f, "pm"];

zl = 2*Diagonal[pm] - 1;
sz = Import[f, "sz"];
(zl - sz // Chop) == Table[0, nSites]
nn = Table[
   1/4 (1 - zl[[i]] - zl[[j]] + zz[[i, j]]), {i, nSites}, {j, 
    nSites}];


(* pm correlations. *)
iCenter = nSites/2;
icShifts = {-45, -30, -15, 0};
dataPM = Table[
   ic = iCenter + i;
   rline = Table[{2*r, pm[[ic - r, ic + r]]}, {r, 2, 15}];
   fit = 
    FindFit[rline, {af + bf*x^(2*z - 1)*(1 + cf*x^z)}, {af, bf, cf}, 
     x];
   eq = af + bf*x^(2*z - 1)*(1 + cf*x^z) /. fit;
   {rline, eq}
   , {i, icShifts}];
{rlines, eqs} = Transpose[dataPM];

(* setup plot. *)
colors = colors0[[;; Length[icShifts]]];
(*leg1=SwatchLegend[
Table[{ic,c}=params;
Style[ic,c],
{params,Transpose[{iCenter+icShifts,colors}]}],
LabelStyle\[Rule]{FontSize\[Rule]14},
LegendFunction\[Rule]frame1,
LegendLayout->{"Column",2},
LegendMarkers->markers];
*)
fig1 = ListPlot[
   rlines,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, FontSize -> 24] & /@ {Abs[i - j], 
       "\[LeftAngleBracket]\!\(\*SubsuperscriptBox[\(S\), \(i\), \
\(+\)]\)\!\(\*SubsuperscriptBox[\(S\), \(j\), \(-\)]\)\
\[RightAngleBracket] \[Times]\!\(\*SuperscriptBox[\(10\), \
\(2\)]\)"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {20, Black, FontFamily -> "Times New Roman"},
   PlotMarkers -> "OpenMarkers",
   PlotRange -> {{1, 32}, {0.009, 0.032}},
   PlotStyle -> Table[{Thickness[0.01], c}, {c, colors}],
   (*PlotLegends->Placed[leg1,{0.82,0.87}],*)
   
   FrameTicks -> {{{{0.010, "1.0"}, {0.015, "1.5"}, {0.020, 
        "2.0"}, {0.025, "2.5"}, {0.030, "3.0"}}, None}, Automatic}
   ];
fig2 = Plot[eqs, {x, 1, 200},
   PlotRange -> All,
   PlotStyle -> Table[{Thickness[0.006], c}, {c, colors}]];
fig3 = Show[{fig1, fig2}]
Export[StringTemplate["figs/pmfit-n``-a``-Jz``-magnon``.pdf"][nSites, 
  alpha, Jz, numMagnons], fig3]


(* nn correlations. *)
iCenter = nSites/2;
icShifts = {-45, -30, -15, 0};
dataZZ = Table[
   ic = iCenter + i;
   rline = Table[{2*r, nn[[ic - r, ic + r]]}, {r, 2, 15}];
   fit = 
    FindFit[rline, {af + bf*x^(2*z - 2)*(1 + cf*x^z)}, {af, bf, cf}, 
     x];
   eq = af + bf*x^(2*z - 2)*(1 + cf*x^z) /. fit;
   {rline, eq}
   , {i, icShifts}];
{rlines, eqs} = Transpose[dataZZ];

(* setup plot. *)
colors = colors0[[;; Length[icShifts]]];
leg1 = SwatchLegend[
   Table[{ic, c} = params;
    Style[ic, c],
    {params, Transpose[{iCenter + icShifts, colors}]}],
   LabelStyle -> {FontSize -> 18},
   LegendFunction -> frame1,
   LegendLayout -> {"Column", 2},
   LegendMarkers -> markers];

fig1 = ListPlot[
   rlines,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, FontSize -> 24] & /@ {Abs[i - j], 
       "\[LeftAngleBracket]\!\(\*SubscriptBox[\(n\), \
\(i\)]\)\!\(\*SubscriptBox[\(n\), \(j\)]\)\[RightAngleBracket] \
\[Times]\!\(\*SuperscriptBox[\(10\), \(4\)]\)"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {20, Black, FontFamily -> "Times New Roman"},
   PlotMarkers -> "OpenMarkers",
   PlotRange -> {{2, 32}, All},
   PlotStyle -> Table[{Thickness[0.01], c}, {c, colors}],
   PlotLegends -> Placed[leg1, {0.78, 0.83}],
   FrameTicks -> {{{{0., "0"}, {0.0002, "2"}, {0.0004, "4"}, {0.0006, 
        "6"}, {0.0008, "8"}, {0.0010, "10"}, {0.0012, "12"}}, None}, 
     Automatic}
   ];
fig2 = Plot[eqs, {x, 2, 200},
   PlotRange -> All,
   PlotStyle -> Table[{Thickness[0.006], c}, {c, colors}]];
fig3 = Show[{fig1, fig2}]
Export[StringTemplate["figs/zzfit-n``-a``-Jz``-magnon``.pdf"][nSites, 
  alpha, Jz, numMagnons], fig3]


(* Part II. Demonstrate universality: Show fitting robustness. *)
(* model parameters. *)
nSites = 200;
alpha = 1.9;
z := alpha - 1;
Jz = 1.9;
(* fitting parameters*)
rLong = 15;
rShortList = {1, 2, 3, 4};
icCut = 40;
nMagnonsList = {2, 3, 4};

(* Universality fitting on different magnon number states. *)
data = Table[
   (* load MPS data. *)
   f = FileNameJoin[{Directory[], "data", 
      StringTemplate["magnon-num_magnons``-n``-a``-Jz``-bc-MPS.h5"][
       numMagnons, nSites, alpha, Jz]}];
   zz = Import[f, "zz"];
   pm = Import[f, "pm"];
   zl = 2*Diagonal[pm] - 1;
   nn = Table[
     1/4 (1 - zl[[i]] - zl[[j]] + zz[[i, j]]), {i, nSites}, {j, 
      nSites}];
   
   (* universal fitting C(X). *)
   rrPM = Table[
     rl = Table[
       rline = 
        Table[{2*r, pm[[ic - r, ic + r]]}, {r, rShort, rLong}];
       fit = 
        FindFit[rline, {af + bf*x^(2*z - 1)*(1 + cf*x^z)}, {af, bf, 
          cf}, x];
       c = (4*Pi*bf)/(Gamma[1 - 2*z]*Sin[Pi*z]) /. fit,
       {rShort, rShortList}];
     {ic, Around[Mean[rl], StandardDeviation[rl]]},
     {ic, 1 + icCut, nSites - icCut, 4}];
   rrNN = Table[
     rl = Table[
       rline = 
        Table[{2*r, nn[[ic - r, ic + r]]}, {r, rShort, rLong}];
       fit = 
        FindFit[rline, {af + bf*x^(2*z - 2)*(1 + cf*x^z)}, {af, bf, 
          cf}, x];
       c = (4*Pi*Pi*bf)/(Gamma[1 - z]^2*Sin[Pi*z/2]^2) /. fit,
       {rShort, rShortList}];
     {ic, Around[Mean[rl], StandardDeviation[rl]]},
     {ic, 3 + icCut, nSites - icCut, 4}];
   
   {rrPM, rrNN}
   , {numMagnons, nMagnonsList}];

(* plot config. *)
colors = colors0[[;; 2*Length[nMagnonsList]]];
colors = Join[colors[[;; 5]], {Gray}];
labels = Join @@ Table[
    {Subscript[
      "\[LeftAngleBracket]\!\(\*SubsuperscriptBox[\(S\), \(i\), \
\(+\)]\)\!\(\*SubsuperscriptBox[\(S\), \(j\), \(-\)]\)\
\[RightAngleBracket]", StringTemplate["``"][numMagnons]], 
     Subscript[
      "\[LeftAngleBracket]\!\(\*SubscriptBox[\(n\), \
\(i\)]\)\!\(\*SubscriptBox[\(n\), \(j\)]\)\[RightAngleBracket]", 
      StringTemplate["``"][numMagnons]]},
    {numMagnons, {2, 3, 4}}];
labelsStyled = Style @@@ Transpose[{labels, colors}];
leg1 = SwatchLegend[
   labelsStyled,
   LabelStyle -> {FontSize -> 12},
   LegendFunction -> frame1,
   LegendLayout -> {"Row", 3},
   (*LegendLayout->{"Column",1},*)
   LegendMarkers -> markers];
aspratio = 4/5;
imgsize = 500;

(* plot universality C(X). *)
fig1 = ListPlot[
   Join @@ data,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, FontSize -> 18] & /@ {X, 
       TraditionalForm[HoldForm[c[X] "/" u^2]]}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {18, Black, FontFamily -> "Times New Roman"},
   (*PlotMarkers->"OpenMarkers",*)
   PlotMarkers -> markers,
   PlotRange -> {{40, 161}, All},
   PlotStyle -> Table[{Thickness[0.006], c}, {c, colors}],
   PlotLegends -> Placed[leg1, {0.83, 0.87}]
   ];
Export[StringTemplate[
   "figs/cxfit_universality-n``-a``-Jz``-magnons234.pdf"][nSites, 
  alpha, Jz], fig1]
Show[fig1]
