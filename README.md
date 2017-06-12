# numerico

Numerics in javascript.

This is an extension of Sebastian Loisels numericjs including improvements from several other sources (which I hope to credit here shortly. If you are uncredited, please reach out). 

Numerico has three main advantages over numericjs:

1. Wraps the numeric library
2. A jacobi algorithm for explicitly diagonalizing Hermitian matrices. The algorithm works problemless even for systems with degenerate eigenvalues.
3. It includes a gulpfile for easy selection and minification of the numeric packages.

It is used in [quantumgraphs.com](https://quantumgraphs.com).
