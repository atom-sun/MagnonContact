On[Assert];

Clear["Global`*"];
Remove["Global`*"];

path0 = DirectoryName[$InputFileName];
path1 = NotebookDirectory[];
If[! DirectoryQ[path0], SetDirectory[path1], SetDirectory[path0]];
Print[Directory[]]


Plot[Sin[x], {x, 0, 10}]


(* Setup plot basics. *)
(*colors=System`PlotThemeDump`$ThemeDefaultGradient*)
Table[{i, ColorData[i] /@ Range[10]}, {i, 115}];
colors0 = ColorData[110] /@ Range[10];
markers = {circle, uptriangle, diamond, square, downtriangle} =
   Charting`CommonDump`GraphicsOpenPlotMarkers[][[;; 5]];
aspratio = 4/5;
imgsize = 400;
frame1[legend_] := Framed[legend,
   Background -> White,
   FrameStyle -> Black,
   FrameMargins -> 0,
   RoundingRadius -> 10
   ];


(* theoretical calculation. *)
alphaList = {1.75, 1.8, 1.85, 1.9};

data = Table[
   z := alpha - 1;
   
   res = Table[
     \[Omega] = i\[Omega];
     roots1 = 
      x /. NSolve[\[Omega] - Abs[1/2 + x]^z - Abs[1/2 - x]^z == 0, x, 
        Reals];
     roots2 = 
      x /. NSolve[\[Omega] - Abs[x]^z - Abs[1 + x]^z == 0, x, Reals];
     
     (*Diagram 1*)
     DtildeInvRe = (-1/2) NIntegrate[
        1/(2 \[Pi]) (1/(\[Omega] - Abs[1/2 + x]^z - Abs[1/2 - x]^z) + 
           1/(2 Abs[x]^z)), {x, -10^6, roots1[[1]], roots1[[2]], 
         10^6}, Method -> {"PrincipalValue"}];
     
     DtildeInvIm = Table[
        If[roots1[[i]] < -1/2, 
         aa = (-1/2) 1/(
            2 \[Pi]) ((-I  \[Pi] )/
             Abs[-z (-(1/2) - x)^(z - 1) - z (1/2 - x)^(z - 1)]) /. 
           x -> roots1[[i]]];
        If[1/2 > roots1[[i]] > -1/2, 
         aa = (-1/2) 1/(
            2 \[Pi]) ((-I  \[Pi] )/
             Abs[z (1/2 + x)^(z - 1) - z (1/2 - x)^(z - 1)]) /. 
           x -> roots1[[i]]];
        If[roots1[[i]] > 1/2, 
         aa = (-1/2) 1/(
            2 \[Pi]) ((-I  \[Pi] )/
             Abs[z (1/2 + x)^(z - 1) + z (x - 1/2)^(z - 1)]) /. 
           x -> roots1[[i]]];
        aa, {i, 1, Length[roots1]}] // Total;
     DtildeInv = DtildeInvRe + DtildeInvIm;
     
     Int1Re = 
      NIntegrate[
       1/(2 \[Pi]) (1/(\[Omega] - Abs[x]^z - Abs[1 + x]^z) 1/(
          2 Abs[x]^z)), {x, -Infinity, roots2[[1]], 0, roots2[[2]], 
        Infinity}, Method -> {"PrincipalValue"}];
     Int1Im = Table[
        If[roots2[[i]] < -1, 
         aa = 1/(2 \[Pi]) ((-I  \[Pi] )/
             Abs[-z (-x)^(z - 1) - z (-1 - x)^(z - 1)]) 1/(
            2 Abs[x]^z) /. x -> roots2[[i]]];
        If[0 > roots2[[i]] > -1, 
         aa = 1/(2 \[Pi]) ((-I  \[Pi] )/
             Abs[-z (-x)^(z - 1) + z (1 + x)^(z - 1)]) 1/(
            2 Abs[x]^z) /. x -> roots2[[i]]];
        If[roots2[[i]] > 0, 
         aa = 1/(2 \[Pi]) ((-I  \[Pi] )/
             Abs[z (x)^(z - 1) + z (1 + x)^(z - 1)]) 1/(
            2 Abs[x]^z) /. x -> roots2[[i]]];
        aa, {i, 1, Length[roots2]}] // Total;
     Int1 = Int1Re + Int1Im;
     I1\[Omega] = -Im[(Int1^2)/DtildeInv ]/\[Pi];
     
     (*Diagram 2*)
     I2\[Omega] = 
      Table[If[roots2[[i]] < -1, 
         aa = 1/(2 \[Pi]) (1/
             Abs[-z (-x)^(z - 1) - z (-1 - x)^(z - 1)]) 1/(
            4 Abs[x]^(2 z)) /. x -> roots2[[i]]];
        If[0 > roots2[[i]] > -1, 
         aa = 1/(2 \[Pi]) (1/
             Abs[-z (-x)^(z - 1) + z (1 + x)^(z - 1)]) 1/(
            4 Abs[x]^(2 z)) /. x -> roots2[[i]]];
        If[roots2[[i]] > 0, 
         aa = 1/(2 \[Pi]) (1/
             Abs[z (x)^(z - 1) + z (1 + x)^(z - 1)]) 1/(
            4 Abs[x]^(2 z)) /. x -> roots2[[i]]];
        aa, {i, 1, Length[roots2]}] // Total;
     
     (*Diagram 3*)
     I3\[Omega] = 
      Table[If[roots2[[i]] < -1, 
         aa = 1/(2 \[Pi]) (1/
             Abs[-z (-x)^(z - 1) - z (-1 - x)^(z - 1)]) 1/(
            4 Abs[x]^z Abs[1 + x]^z) /. x -> roots2[[i]]];
        If[0 > roots2[[i]] > -1, 
         aa = 1/(2 \[Pi]) (1/
             Abs[-z (-x)^(z - 1) + z (1 + x)^(z - 1)]) 1/(
            4 Abs[x]^z Abs[1 + x]^z) /. x -> roots2[[i]]];
        If[roots2[[i]] > 0, 
         aa = 1/(2 \[Pi]) (1/
             Abs[z (x)^(z - 1) + z (1 + x)^(z - 1)]) 1/(
            4 Abs[x]^z Abs[1 + x]^z) /. x -> roots2[[i]]];
        aa, {i, 1, Length[roots2]}] // Total;
     
     I\[Omega] = I1\[Omega] + I2\[Omega] + I3\[Omega];
     {\[Omega], I\[Omega]}, {i\[Omega], 1/2^(z - 1) + 0.01, 2, 0.01}]
   
   , {alpha, alphaList}];


(* plot I(omega) vs omega. *)
colors = colors0[[;; Length[alphaList]]];
labels = Table[
   {alpha, c} = params;
   Style[StringTemplate["\[Alpha]=``"][alpha], c],
   {params, Transpose[{alphaList, colors}]}];
leg1 = SwatchLegend[
   labels,
   LabelStyle -> {FontSize -> 12},
   LegendFunction -> frame1,
   LegendLayout -> {"Column", 1},
   LegendMarkers -> markers];

fig1 = ListLinePlot[
   data,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, FontSize -> 18] & /@ {TraditionalForm[
        HoldForm[\[Omega] "/" u*k^z]], 
       TraditionalForm[
        HoldForm[u^3 k^(3*z - 1) \[ScriptCapitalI][\[Omega], k]]]}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {18, Black, FontFamily -> "Times New Roman"},
   PlotMarkers -> "OpenMarkers",
   PlotRange -> {All, {-0.01, 0.6}},
   PlotStyle -> Table[{Thickness[0.006], c}, {c, colors}],
   PlotLegends -> Placed[leg1, {0.86, 0.8}]
   ];
Export["figs/iomega.pdf", fig1];
Show[fig1]
