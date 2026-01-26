To compile:
gfortran -c -w rmatrix.f trb8.f

gfortran hormat-so.f rmatrix.o trb8.o -o hormat-so -llapack -lblas 

To run:
./hormat-so <data-file

I added some testing loop in hormat around the wfn writing section and tested it with 2fnr.
2fnr can read it now (after I take out the duplicated L=0 one) but the resulting (p,t) xsec is very wrong.

I just need to compare the wfns from Greg's output to the ones from hormat.
Natasha said she found an error in Greg's code for the TPM, the imaginary radius. Need to look into.