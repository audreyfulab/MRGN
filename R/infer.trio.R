


#' A function to infer the causal network of a trio
#'
#' This function takes in a matrix representing a trio (with or without confounding variables) and performs the
#' causal network inference. It wraps the functions get.freq(), Reg(), PermReg(), and class.vec() to infer an input trio.
#' It returns a one dimensional dataframe containing the results of the coefficient and marginal tests and the
#' inferred model structure.
#' If the genetic variant is copy number alteration (CNA), permutation tests may 
#' be appropriate. To perform such tests, permutation tests, set is.CNA=TRUE and use.perm=TRUE and specify nperms.
#'
#' @param trio A dataframe with at least 3 columns. The first column should be the genetic variant, and the second
#'             and third columns the molecular phenotypes. Data in columns 4+ are treated as confounding variables
#' @param use.perm (logical) if TRUE the permutation test is used for trios with minor allele freq < gamma. If FALSE
#'                 test is ignored for all trios (default = FALSE)
#' @param gamma The minor allele frequency threshold for which the permutation test should be used such that permutation is
#'              performed when minor allele freq < gamma (when use.perm = TRUE).
#' @param is.CNA (logical) for use when use.perm=TRUE and the genetic variant is a copy number alteration so that the
#'               threshold gamma is ignored. default = FALSE
#' @param alpha The rejection threshold for all Wald tests (default = 0.01)
#' @param nperms The number of permutations within each genotype (default = 1,000)
#' @param compute.nominal (logical) Set to FALSE, so that ermutation p-values are computed as
#'                        (B+1)/(nperms+1) where B is the number of permuted statistics >= the observed. A peudocount of 1 is added to the calculation.
#' @param verbose (logical) if TRUE results of the regressions are printed
#' @examples
#' #inference on a single trio with an eQTL as genetic variant
#' result=infer.trio(M1trio)
#' print(result)
#'
#' \dontrun{
#' #fast example on 10 eQTL trios from the built in dataset WBtrios
#' results = sapply(WBtrios[1:10], function(x) infer.trio(x))
#' print(results)
#' #return just the inferred model topology
#' models = sapply(WBtrios[1:10], function(x) infer.trio(x)$Inferred.Model)
#' print(models)
#' 
#' #inference on a single trio with a copy number alteration (CNA)
#' #as genetic variant
#' dim (CNAtrios[[1]])
#' # 564x37
#' result=infer.trio(CNAtrios[1], is.CNA = TRUE, use.perm = TRUE)
#' print(result)
#' #fast example on 10 CNA trios from the built in dataset CNAtrios using permutation
#' models = sapply(CNAtrios[1:10], function(x) infer.trio(x, is.CNA = TRUE, use.perm = TRUE, nperms = 1000)$Inferred.Model)
#' print(models)
#' }
#' @return a dataframe of dimension 1 x 18 with the following columns:
#'   \describe{
#'   \item{b11}{binary value for the conditional test T1 ~ V | T2,U. 1 for rejection of the null hypothesis, and 0 for no rejection.}
#'   \item{b12}{binary value for the conditional test T1 ~ T2 | V,U}
#'   \item{b21}{binary value for the conditional test T2 ~ V | T1,U}
#'   \item{b22}{binary value for the conditional test T2 ~ T1 | V,U}
#'   \item{V1:T1}{binary value for the marginal test between V1 and T1}
#'   \item{V1:T2}{binary value for the marginal test between V1 and T2}
#'   \item{pb11}{the p-value for the conditional test T1 ~ V | T2,U}
#'   \item{pb12}{the p-value for the conditional test T1 ~ T2 | V,U}
#'   \item{pb21}{the p-value for the conditional test T2 ~ V | T1,U}
#'   \item{pb22}{the p-value for the conditional test T2 ~ T1 | V,U}
#'   \item{pV1:T1}{the p-value for the marginal test between V1 and T1}
#'   \item{pV1:T2}{the p-value for the marginal test between V1 and T2}
#'   \item{Minor.freq}{the calculated frequency of the minor allele of the genetic variant. Ignore when 
#'  the genetic variant is CNA}
#'   \item{coef11}{the coefficient estimate of beta_11}
#'   \item{coef12}{the coefficient estimate of beta_12}
#'   \item{coef21}{the coefficient estimate of beta_21}
#'   \item{coef22}{the coefficient estimate of beta_22}
#'   \item{Inferred.Model}{a string indicating the inferred model type as returned by \eqn{class.vec()}}
#'   }
#'
#' @export infer.trio



####################################################################
#a wrapper function for get.freq(), Reg(), and PermReg() to infer the trio
#combines the functions from sections 1.1-1.2
infer.trio=function(trio=NULL, use.perm = FALSE, gamma=0.05, is.CNA = FALSE, alpha=0.01, nperms=1000, compute.nominal=FALSE, verbose=FALSE){

  #ensure trio is a dataframe for later functions:
  trio = as.data.frame(trio)
  #preallocate indicator vectors
  xp=NULL
  rp=NULL

  #preform the standard regressions and outputs t-stat and p-values
  #input is a trio with the variant in the first column followed by genes and confounders [V,T1,T2,U]
  #step 1
  pt.out=Reg(data = trio, verbose=verbose)



  #check the frequency of the minor allele using get.freq()
  #step 2
  minor=get.freq(V=trio[,1], is.CNA = is.CNA)

  #step 2.1
  if(use.perm == TRUE & (minor<gamma | is.CNA)){
    #preform permuted regression (section 1.2) for rare variants
    #note: when compute.nominal=FALSE, permutation p-values use the (B+1)/(m+1) correction
    #(Phipson & Smyth, 2010) to include the observed statistic in the reference distribution
    pvals=PermReg(trio = na.omit(trio),
                  t.obs12 = pt.out$tvals[2],
                  t.obs22 = pt.out$tvals[4],
                  p11 = pt.out$pvals[1],
                  p21 = pt.out$pvals[3],
                  m = nperms,
                  compute.nominal = compute.nominal)

  }else{
    #else return the pvals from standard reg.
    pvals=pt.out$pvals

  }

  #---steps 3-4
  #convert xp to indicator vector
  #section 1.1 step
  xp=ifelse(pvals<alpha, 1, 0)

  #preform marginal tests
  cors=c(stats::cor.test(trio[,1], trio[,2],use="pairwise.complete.obs")$p.value,
         stats::cor.test(trio[,1], trio[,3],use="pairwise.complete.obs")$p.value)
  #convert to indicator vector
  rp=ifelse(cors<alpha, 1, 0)
  #combine all useful stats - add indicator
  all.stats=c(append(xp, rp), append(pvals, cors), minor, pt.out$tvals)

  names(all.stats)=c("b11","b12", "b21","b22", "V1:T1", "V1:T2", "pb11",
                     "pb12", "pb21","pb22","pV1:T1","pV1:T2", "Minor.freq", 
                     "coef11", "coef12", "coef21", "coef22")
  all.stats = as.data.frame(t(all.stats))
  all.stats$Inferred.Model = MRGN::class.vec(all.stats)

  return(all.stats)
}
