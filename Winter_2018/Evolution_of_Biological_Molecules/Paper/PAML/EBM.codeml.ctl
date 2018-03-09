seqfile = /home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/Paper/FLS2_Genes_Multisp_aligned.fst
outfile = /home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/Paper/PAML/temp.out
treefile = /home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/Paper/bootstrap.tree
verbose = 2
noisy = 9
runmode = 0 * Use user-supplied tree? 0:yes 1:no
seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
aaRatefile = /home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/Paper/PAML/jones.dat * only used for aa seqs with model=empirical(_F) dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
model = 2 * 2 = model-given frequencies 3 = model+estimated_frequencies
fix_alpha = 0 * 0: estimate gamma shape parameter; 1: fix it at alpha alpha = 0. initial or fixed alpha, 0:infinity (constant rate)
alpha = 0.5 * initial or fixed alpha (fixed if fix_alpha = 1; intitial alpha if fix_alpha = 0)
RateAncestor = 1 * infer rates or ancestor? (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
cleandata = 0 * remove sites with alignment gaps (1:yes, 0:no)? 
fix_blength = 0 * user-supplied branch lengths 0: ignore, -1: random, 1: initial, 2: fixed
method = 1   * branch-swapping 0: simultaneous; 1: one branch at a time
aaDist = 0 * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a 7:AAClasses
Small_Diff = .5e-6
Malpha = 0  * different alphas for genes
ncatG = 3  * # of categories in dG of NSsites models
fix_rho = 1  * 0: estimate rho; 1: fix it at rho
rho = 0. * initial or fixed rho,   0:no correlation
getSE = 0 * 0: don't want them, 1: want S.E.s of estimates 
