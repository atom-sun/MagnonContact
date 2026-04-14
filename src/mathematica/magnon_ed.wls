On[Assert];

Clear["Global`*"];
Remove["Global`*"];

path0 = DirectoryName[$InputFileName];
path1 = NotebookDirectory[];
If[! DirectoryQ[path0], SetDirectory[path1], SetDirectory[path0]];
Print[Directory[]]


(* ED: one-magnon subspace. *)
(* define Hamiltonian terms in one-magnon subspaces. *)
hop1[nSites_Integer, alpha_Real] := 
  Table[If[l == m, 0, -2./Abs[l - m]^alpha], {l, nSites}, {m, nSites}];
hzz1[nSites_Integer, boundaryCompensation_Symbol : True] := 
  DiagonalMatrix[
   Join[{nSites - 3 - 2*Boole[boundaryCompensation]}, 
    Table[nSites - 5, 
     nSites - 2], {nSites - 3 - 2*Boole[boundaryCompensation]}]];

(* model parameters. *)
nSites = 20;
alpha = 1.75;
Jz = 2.5;

(* solving *)
{e0, v0, z, zz, pm} = Module[
   {ham, zjOp, zizjOp, SpmOp, erg, vec, where, e0, v0, z, zz, pm},
   
   (* model Hamiltonian. *)
   ham = hop1[nSites, alpha] + (-Jz)*hzz1[nSites, True];
   ham // MatrixForm;
   
   (* Correlation functinos. *)
   zjOp[j_Integer] :=
    DiagonalMatrix[Table[1 - 2*Boole[l == j], {l, nSites}]];
   zizjOp[i_Integer, j_Integer] /; i != j :=
    DiagonalMatrix[
     Table[1 - 2*Boole[l == i] - 2*Boole[l == j], {l, nSites}]];
   zizjOp[i_Integer, j_Integer] /; i == j := IdentityMatrix[nSites];
   SpmOp[i_Integer, j_Integer] /; i != j :=
    Table[Boole[l == j]*Boole[m == i], {l, nSites}, {m, nSites}];
   SpmOp[i_Integer, j_Integer] /; i == j :=
    DiagonalMatrix[Table[1 - Boole[l == i], {l, nSites}]];
   
   (* solving. *)
   {erg, vec} = Eigensystem[N[ham]];
   where = 
    DeleteDuplicates[
     Flatten[Map[Position[erg, #] &, Sort[erg, Less]]]];
   erg = erg[[where]];
   vec = vec[[where]];
   
   (* ground state. *)
   e0 = erg[[1]];
   v0 = vec[[1]];
   Print[e0];
   (*ListPlot[v0];*)
   
   (* compute expects. *)
   z = Table[vec[[1]] . zjOp[j] . vec[[1]], {j, nSites}];
   Print[Total[z]];
   zz = Table[v0 . zizjOp[i, j] . v0, {i, nSites}, {j, nSites}];
   pm = Table[v0 . SpmOp[i, j] . v0, {i, nSites}, {j, nSites}];
   
   {e0, v0, z, zz, pm}
   ];

Print[e0]
ListPlot[v0]
Print[z // Total]
zz[[;; 5, ;; 5]] // MatrixForm
pm[[;; 5, ;; 5]] // MatrixForm


(* ED: two-magnons subspace. *)
(* illustration: states arrangement. *)
nSites = 5;
pig = Table[
   Subscript[a, k, l, m, n],
   {k, nSites - 1}, {l, k + 1, nSites}, {m, nSites - 1}, {n, m + 1, 
    nSites}];
Join @@@ (Join @@ pig) // MatrixForm

(* define Hamiltonian terms in two-magnon subspaces. *)
hop2[nSites_Integer, alpha_Real] :=
  Join @@@ (Join @@ Table[
      Which[
       k == l, -2./Abs[m - j]^alpha,
       j == m, -2./Abs[k - l]^alpha,
       j == l && k != m, -2./Abs[k - m]^alpha,
       j != l && k == m, -2./Abs[j - l]^alpha,
       True, 0],
      {j, nSites - 1}, {k, j + 1, nSites}, {l, nSites - 1}, {m, l + 1,
        nSites}]);
hzz2[nSites_Integer, boundaryCompensation_Symbol : True] :=
  DiagonalMatrix[
   (nSites - 1) - 2*Join @@ Table[
       Boole[l > 1] + 2*Boole[m - l > 1] + Boole[m < nSites],
       {l, nSites - 1}, {m, l + 1, nSites}]
    + Boole[boundaryCompensation]*(-2)*
     Join @@ Table[
       Boole[l == 1] + Boole[m == nSites], {l, nSites - 1}, {m, l + 1,
         nSites}]
   ];

(* model parameters. *)
nSites = 20;
alpha = 1.75;
Jz = 2.5;

(* solving. *)
{e0, v0, z, zz, pm} = Module[
   {nStates, ham, zjOp, zizjOp, SpmOp, erg, vec, where, e0, v0, z, zz, pm},
   
   nStates = nSites*(nSites - 1)/2;
   
   ham = hop2[nSites, alpha] + (-Jz)*hzz2[nSites, True];
   ham // MatrixForm;
   
   (* Correlation functinos. *)
   zjOp[j_Integer] :=
    DiagonalMatrix[Join @@ Table[1 - 2*Boole[j == l || j == m],
       {l, nSites - 1}, {m, l + 1, nSites}]];
   zizjOp[i_Integer, j_Integer] /; i != j :=
    DiagonalMatrix[Join @@ Table[
       1 - 
        2*Mod[Boole[l == i] + Boole[l == j] + Boole[m == i] + 
           Boole[m == j], 2],
       {l, nSites - 1}, {m, l + 1, nSites}]];
   zizjOp[i_Integer, j_Integer] /; i == j := IdentityMatrix[nStates];
   SpmOp[i_Integer, j_Integer] /; i != j :=
    Join @@@ (Join @@ Table[
        Boole[
         (k == j && l != i && m == i && n != j && l == n) ||
          (k == j && l != i && n == i && m != j && l == m) ||
          (l == j && k != i && m == i && n != j && k == n) ||
          (l == j && k != i && n == i && m != j && k == m)
         ],
        {k, nSites - 1}, {l, k + 1, nSites}, {m, nSites - 1}, {n, 
         m + 1, nSites}]);
   SpmOp[i_Integer, j_Integer] /; i == j :=
    DiagonalMatrix[Join @@ Table[1 - Boole[l == i || m == i],
       {l, nSites - 1}, {m, l + 1, nSites}]];
   
   (* solving. *)
   {erg, vec} = Eigensystem[N[ham]];
   where = 
    DeleteDuplicates[
     Flatten[Map[Position[erg, #] &, Sort[erg, Less]]]];
   erg = erg[[where]];
   vec = vec[[where]];
   
   (* ground state. *)
   e0 = erg[[1]];
   v0 = vec[[1]];
   Print[e0];
   (*ListPlot[v0];*)
   
   (* compute expects. *)
   z = Table[vec[[1]] . zjOp[j] . vec[[1]], {j, nSites}];
   Print[Total[z]];
   zz = Table[v0 . zizjOp[i, j] . v0, {i, nSites}, {j, nSites}];
   pm = Table[v0 . SpmOp[i, j] . v0, {i, nSites}, {j, nSites}];
   
   {e0, v0, z, zz, pm}
   ];

Print[e0]
ListPlot[v0]
Print[z // Total]
zz[[;; 5, ;; 5]] // MatrixForm
pm[[;; 5, ;; 5]] // MatrixForm
