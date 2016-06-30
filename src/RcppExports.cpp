// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// DistKd
void DistKd(SEXP Rx, SEXP Ry, SEXP RPointWeight, SEXP RWeight, SEXP RDist, SEXP RIsReferenceType, SEXP RIsNeighborType);
RcppExport SEXP dbmss_DistKd(SEXP RxSEXP, SEXP RySEXP, SEXP RPointWeightSEXP, SEXP RWeightSEXP, SEXP RDistSEXP, SEXP RIsReferenceTypeSEXP, SEXP RIsNeighborTypeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type Rx(RxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Ry(RySEXP);
    Rcpp::traits::input_parameter< SEXP >::type RPointWeight(RPointWeightSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RWeight(RWeightSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RDist(RDistSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RIsReferenceType(RIsReferenceTypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RIsNeighborType(RIsNeighborTypeSEXP);
    DistKd(Rx, Ry, RPointWeight, RWeight, RDist, RIsReferenceType, RIsNeighborType);
    return R_NilValue;
END_RCPP
}
// CountNbdKd
void CountNbdKd(SEXP Rr, SEXP Rx, SEXP Ry, SEXP RWeight, SEXP RNbd, SEXP RIsReferenceType, SEXP RIsNeighborType);
RcppExport SEXP dbmss_CountNbdKd(SEXP RrSEXP, SEXP RxSEXP, SEXP RySEXP, SEXP RWeightSEXP, SEXP RNbdSEXP, SEXP RIsReferenceTypeSEXP, SEXP RIsNeighborTypeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type Rr(RrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rx(RxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Ry(RySEXP);
    Rcpp::traits::input_parameter< SEXP >::type RWeight(RWeightSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNbd(RNbdSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RIsReferenceType(RIsReferenceTypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RIsNeighborType(RIsNeighborTypeSEXP);
    CountNbdKd(Rr, Rx, Ry, RWeight, RNbd, RIsReferenceType, RIsNeighborType);
    return R_NilValue;
END_RCPP
}
// parallelCountNbd
NumericMatrix parallelCountNbd(NumericVector r, NumericVector x, NumericVector y, NumericVector Weight, LogicalVector IsReferenceType, LogicalVector IsNeighborType);
RcppExport SEXP dbmss_parallelCountNbd(SEXP rSEXP, SEXP xSEXP, SEXP ySEXP, SEXP WeightSEXP, SEXP IsReferenceTypeSEXP, SEXP IsNeighborTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Weight(WeightSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsReferenceType(IsReferenceTypeSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsNeighborType(IsNeighborTypeSEXP);
    __result = Rcpp::wrap(parallelCountNbd(r, x, y, Weight, IsReferenceType, IsNeighborType));
    return __result;
END_RCPP
}
// parallelCountNbdDt
NumericMatrix parallelCountNbdDt(NumericVector r, NumericMatrix Dmatrix, NumericVector Weight, LogicalVector IsReferenceType, LogicalVector IsNeighborType);
RcppExport SEXP dbmss_parallelCountNbdDt(SEXP rSEXP, SEXP DmatrixSEXP, SEXP WeightSEXP, SEXP IsReferenceTypeSEXP, SEXP IsNeighborTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dmatrix(DmatrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Weight(WeightSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsReferenceType(IsReferenceTypeSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsNeighborType(IsNeighborTypeSEXP);
    __result = Rcpp::wrap(parallelCountNbdDt(r, Dmatrix, Weight, IsReferenceType, IsNeighborType));
    return __result;
END_RCPP
}
// parallelCountNbdCC
NumericMatrix parallelCountNbdCC(NumericVector r, NumericVector x, NumericVector y, NumericVector Weight, LogicalVector IsReferenceType, LogicalVector IsNeighborType);
RcppExport SEXP dbmss_parallelCountNbdCC(SEXP rSEXP, SEXP xSEXP, SEXP ySEXP, SEXP WeightSEXP, SEXP IsReferenceTypeSEXP, SEXP IsNeighborTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Weight(WeightSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsReferenceType(IsReferenceTypeSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsNeighborType(IsNeighborTypeSEXP);
    __result = Rcpp::wrap(parallelCountNbdCC(r, x, y, Weight, IsReferenceType, IsNeighborType));
    return __result;
END_RCPP
}
// parallelCountNbdDtCC
NumericMatrix parallelCountNbdDtCC(NumericVector r, NumericMatrix Dmatrix, NumericVector Weight, LogicalVector IsReferenceType, LogicalVector IsNeighborType);
RcppExport SEXP dbmss_parallelCountNbdDtCC(SEXP rSEXP, SEXP DmatrixSEXP, SEXP WeightSEXP, SEXP IsReferenceTypeSEXP, SEXP IsNeighborTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dmatrix(DmatrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Weight(WeightSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsReferenceType(IsReferenceTypeSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type IsNeighborType(IsNeighborTypeSEXP);
    __result = Rcpp::wrap(parallelCountNbdDtCC(r, Dmatrix, Weight, IsReferenceType, IsNeighborType));
    return __result;
END_RCPP
}
// parallelCountNbdm
NumericMatrix parallelCountNbdm(NumericVector x, NumericVector y, IntegerVector ReferencePoints);
RcppExport SEXP dbmss_parallelCountNbdm(SEXP xSEXP, SEXP ySEXP, SEXP ReferencePointsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ReferencePoints(ReferencePointsSEXP);
    __result = Rcpp::wrap(parallelCountNbdm(x, y, ReferencePoints));
    return __result;
END_RCPP
}
