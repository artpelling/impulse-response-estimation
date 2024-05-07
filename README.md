# Compilation
The LaTeX document can be compiled continuously with

``` shell
make latex
```
which will create a symbolic link at `build.pdf` to the latest build.
To only build the document once, execute

``` shell
make document
```
to create `document.pdf`.

# Plots
Any required plots should be generated in `code/plots.py`. The code gets executed and compiled when building with the Makefile. All plots will be saved in `.pdf` format under the [/latex/figures](latex/figures) directory.

# LaTeX style
It is preferred to use notational macros for all mathematical objects s.t. the notation can be easily changed at a later point. The notation is defined in [mynotation.sty](latex/mynotation.sty). Further, [mymathmacros.sty](latex/mymathmacros.sty) contains the setup of all math-related packages as well as optional convenience macros like `\diag{1,2,3}` or set definitions `\R[m][n]` equals `\mathbb{R}^{m\times n}`. Note that `\i` and `\d` remove the math environment and should be used as imaginary unit and in integrals, respectively.

Also, `\( ... \)` is preferred over `$ ... $` for inline math environments.
