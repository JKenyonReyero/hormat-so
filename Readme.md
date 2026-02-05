To compile:
gfortran -c -w rmatrix.f trb8.f

gfortran hormat-so.f rmatrix.o trb8.o -o hormat-so -llapack -lblas 

To run:
./hormat-so <data-file

Output files:

fort.101 - Unknown/not written here yet
fort.102 - Unknown/not written here yet
fort.103 - Unknown/not written here yet
fort.105 - Unknown/not written here yet

fort.111 - Wave functions for LELO
fort.112 - Wave functions for LENLO
fort.113 - Wave functions for ENL

fort.123 - Wave function (just L=0) for ENL

fort.124 - elastic cross section. Use xmgrace -nxy fort.124 to plot it
    1st is Local Equivalent Leading Order (Black)
    2nd is Local Equivalent Next to Leading Order (Red)
    3rd is Exact Non-Local (Green)

fort.211 - Wave functions for LELO   - For use in TwoFNR reading
fort.212 - Wave functions for LENLO  - For use in TwoFNR reading
fort.213 - Wave functions for ENL    - For use in TwoFNR reading

fort.311 - Abs of wave functions for LELO
fort.312 - Abs of wave functions for LENLO
fort.313 - Abs of wave functions for ENL

fort.323 - Abs of wave function (just L=0) for ENL

fort.511 - Unknown/not written here yet
fort.512 - Unknown/not written here yet
fort.802 - Unknown/not written here yet
fort.805 - Unknown/not written here yet

I added some testing loop in hormat around the wfn writing section and tested it with 2fnr.
2fnr can read it now (after I take out the duplicated L=0 one) but the resulting (p,t) xsec is very wrong.

I just need to compare the wfns from Greg's output to the ones from hormat.
Natasha said she found an error in Greg's code for the TPM, the imaginary radius. Need to look into.

28/01/2026
Natasha fixed the wave function writting, so those modifications are included.
Also added more writting output files, written above.
Going to start making (p,t) x secs with twoFNR
If I make the Vso = -0.0001 it will make the twoFNR wfn writting work properly.
Very dirty and quick change but it should work.