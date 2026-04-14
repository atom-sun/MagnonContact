On[Assert];

Clear["Global`*"];
Remove["Global`*"];

path0 = DirectoryName[$InputFileName];
path1 = NotebookDirectory[];
If[! DirectoryQ[path0], SetDirectory[path1], SetDirectory[path0]];
Print[Directory[]]


(* model parameters. *)
nSites = 200;
alpha = 1.75;
Jz = 2.45;

(* numMagnons: 2, 3, 4, 5. *)
numMagnons = 2;

f = FileNameJoin[{Directory[], "data", 
   StringTemplate["magnon-num_magnons``-n``-a``-Jz``-bc-MPS.h5"][
    numMagnons, nSites, alpha, Jz]}]

zz = Import[f, "zz"];
pm = Import[f, "pm"];
