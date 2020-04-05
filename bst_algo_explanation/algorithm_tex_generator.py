with open('algo.tex', 'wt') as f:
    f.write("""
\\documentclass{article}
\\usepackage[utf8]{inputenc}
\\usepackage{pgf}
\\usepackage{tikz}
\\usetikzlibrary{arrows,automata,positioning}
\\usepackage{polski}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{amsfonts}
\\begin{document}

\\begin{figure}[ht]
\centering
\\begin{tikzpicture}[y=0.80pt, x=0.80pt, yscale=1, xscale=1, inner sep=0pt, outer sep=0pt, 
->,>=stealth',shorten >=1pt,auto, thick, trans/.style={thick,->,>=stealth}]
  \\node[state, draw=none] (X)                    {};
""")

    for i in range(10):
        for j in range(10):
            if j == 0:
                if i == 0:
                    f.write(f'\\node[state](I{i}{j})[below = 20 of X] {{{i}, {j}}};\n')
                else:
                    f.write(f'\\node[state](I{i}{j})[below = 20 of I{i - 1}{j}] {{{i}, {j}}};\n')
            else:
                f.write(f'\\node[state](I{i}{j})[right = 20 of I{i}{j - 1}] {{{i}, {j}}};\n')

    f.write("""
    \\path 
(I11) edge [loop above] node {} (I11)
(I11) edge [loop above] node {} (I11)
;
\\end{tikzpicture}
\\label{fig:figure2}
\\end{figure}
\\end{document}
    """)