facts = {1:2,3:4}; % cell array con i factors per ogni strato, in questo
   % caso split-split-plot con un factor nel primo, due nel secondo e quattro
   % nel terzo (per un blocked design, lasci lo strato vuoto, es. {[],1:3})
units = [12,4]; % numero di unita' per ogni strato, in totale 2 x 4 x 3
   % esperimenti
levels = 3; % numero di livelli per ogni variabile, puo' anche essere un
   % vettore con un elemento per ogni variabile se vuoi un numero di
   % livelli diverso per ogni variabile
etas = [1,1]; % rapporto delle varianze fra due strati consecutivi, sono
   % sempre #strati - 1
% criteria = {'D','Nt', 'Pwp', 'Psp'};
model = 'quadratic'; % 'quadratic' 'main' o 'interaction'


example = [
1	-1	1	0
1	-1	-1	1
1	-1	-1	1
1	-1	-1	1
1	-1	1	0
1	-1	-1	1
1	-1	-1	1
1	-1	-1	1
1	1	1	1
1	1	0	1
1	1	-1	0
1	1	1	-1
1	1	1	1
1	1	0	1
1	1	-1	0
1	1	1	-1
-1	-1	0	1
-1	-1	1	-1
-1	-1	1	-1
-1	-1	0	1
0	-1	0	-1
0	-1	1	1
0	-1	-1	-1
0	-1	-1	0
1	1	0	1
1	1	1	-1
1	1	1	1
1	1	-1	0
-1	1	-1	1
-1	1	-1	-1
-1	1	1	-1
-1	1	1	1
0	-1	0	-1
0	-1	1	1
0	-1	-1	0
0	-1	-1	-1
0	-1	0	-1
0	-1	1	1
0	-1	-1	-1
0	-1	-1	0
0	0	1	0
0	0	1	0
0	0	1	0
0	0	1	0
-1	1	1	-1
-1	1	1	1
-1	1	-1	-1
-1	1	-1	1
];
%criteria = {'Psp'};
criteria = {'D','Nt', 'Pwp', 'Psp'};
msopt = MSOpt(facts,units,levels,etas,criteria,model);
%
msopt.Score(example)
