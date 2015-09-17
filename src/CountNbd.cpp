#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void CountNbd(SEXP Rr, SEXP Rx, SEXP Ry, SEXP RWeight, SEXP RNbd, SEXP RIsReferenceType, SEXP RIsNeighborType) {
// Count the number of neighbors around each point
// Only ReferenceType points are considered. 
// The weights of NeighborType points and all points around them are counted.

  //Distances
  NumericVector r(Rr);
  // x, y coordinates of points
  NumericVector x(Rx);
  NumericVector y(Ry);
  // Point weights
  NumericVector Weight(RWeight);
  // Matrix counting the number of neighbors. Modified by this routine.
  NumericMatrix Nbd(RNbd);
  // Boolean vectors describing reference and neighbor points
  IntegerVector IsReferenceType(RIsReferenceType);
  IntegerVector IsNeighborType(RIsNeighborType);

  double Distance2;
  double Nr = r.length();
  NumericVector r2 = r*r;

  for (int i=0; i < (Nbd.nrow()-1); i++) {
    // Consider reference type points
    if (IsReferenceType[i]) {
      // Point j is a neighbor of i. No neighbor is ignored.
      for (int j=i+1; j < Nbd.nrow(); j++) {
        // Calculate squared distance
        Distance2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
        // Ignore point j if it is too far from point i
        if (Distance2 <= r2[Nr-1]) {
          // Find the column of the matrix corresponding to the distance
          int k = 0; 
          while (Distance2 > r2[k]) {
            k++;
          }
          // Add j's weight to i's neighborhood
          Nbd(i, Nr+k) += Weight[j];
          // The neighbor is a point of interest
          if (IsNeighborType[j]) {
            Nbd(i, k) += Weight[j];
          }
          // j is a reference point: add i's neighborhood to it
          if (IsReferenceType[j]) {
            Nbd(j, Nr+k) += Weight[i];
            // i is a point of interest around j
            if (IsNeighborType[i]) {
              Nbd(j, k) += Weight[i];
            }
          }        
        }
      }
    } else {
      // Point i is not a reference point
      for (int j=i+1; j < Nbd.nrow(); j++) {
        // If point j is a reference point, it may be in its neighborhood
        if (IsReferenceType[j]) {
          // Calculate squared distance
          Distance2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
          // Ignore point j if it is too far from point i
          if (Distance2 <= r2[Nr-1]) {
            // Find the column of the matrix corresponding to the distance
            int k = 0; 
            while (Distance2 > r2[k]) {
              k++;
            }
            // Add i's weight to j's neighborhood
            Nbd(j, Nr+k) += Weight[i];
            // i is a point of interest around j
            if (IsNeighborType[i]) {
              Nbd(j, k) += Weight[i];
            }
          }
        }
      }
    }
  }
}

// [[Rcpp::export]]
void CountNbdCC(SEXP Rr, SEXP Rx, SEXP Ry, SEXP RWeight, SEXP RNbd, SEXP RIsReferenceType, SEXP RIsNeighborType) {

  //Distances
  NumericVector r(Rr);
  // x, y coordinates of points
  NumericVector x(Rx);
  NumericVector y(Ry);
  // Point weights
  NumericVector Weight(RWeight);
  // Matrix counting the number of neighbors. Modified by this routine.
  NumericMatrix Nbd(RNbd);
  // Boolean vectors describing reference and neighbor points
  IntegerVector IsReferenceType(RIsReferenceType);
  IntegerVector IsNeighborType(RIsNeighborType);

  double Distance2;
  double Nr = r.length();
  NumericVector r2 = r*r;

  for (int i=0; i < (Nbd.nrow()-1); i++) {
    // Consider cases
    if (IsReferenceType[i]) {
      // Point j is a neighbor of i
      for (int j=i+1; j < Nbd.nrow(); j++) {
        // Ignore point j if it is neither a case nor a control
        if (IsNeighborType[j] || IsReferenceType[j]) {
          // Calculate squared distance
          Distance2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
          // Ignore point j if it is too far from point i
          if (Distance2 <= r2[Nr-1]) {
            // Find the column of the matrix corresponding to the distance
            int k = 0; 
            while (Distance2 > r2[k]) {
              k++;
            }
            // The neighbor is a control: Add j's weight to i's neighborhood
            if (IsNeighborType[j]) {
              Nbd(i, Nr+k) += Weight[j];
            }
            // The neighbor is a case: add it to j's neighborhood (and symetrically add i to j's neighborhood)
            if (IsReferenceType[j]) {
              Nbd(i, k) += Weight[j];
              Nbd(j, k) += Weight[i];
            }
          }
        }
      }
    }
    // Consider controls
    if (IsNeighborType[i]) {
      // Point j is a neighbor of i
      for (int j=i+1; j < Nbd.nrow(); j++) {
        // The neighbor is a case: it is retained because i is a control in the circle around j
        if (IsReferenceType[j]) {
          // Calculate squared distance
          Distance2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
          // Ignore point j if it is too far from point i
          if (Distance2 <= r2[Nr-1]) {
            // Find the column of the matrix corresponding to the distance
            int k = 0; 
            while (Distance2 > r2[k]) {
              k++;
            }
            // Add i's weight to j's neighborhood
            Nbd(j, k) += Weight[i];
          }
        }
      }  
    }
  }
}

// [[Rcpp::export]]
void DistKd(SEXP Rx, SEXP Ry, SEXP RPointWeight, SEXP RWeights, SEXP RDist, SEXP RIsReferenceType, SEXP RIsNeighborType) {
// Fill a distance vector between each (reference point, neighbor point) pair for Kd
// If weighted, also fill a similar vector with the product of point weights

  // x, y coordinates of points
  NumericVector x(Rx);
  NumericVector y(Ry);
  // Point weights
  NumericVector PointWeight(RPointWeight);
  // A vector for weights of point pairs
  NumericVector Weights(RWeights);
  // A vector for distances between point pairs
  NumericVector Dist(RDist);
  // Boolean vectors describing reference and neighbor points
  IntegerVector IsReferenceType(RIsReferenceType);
  IntegerVector IsNeighborType(RIsNeighborType);
  
  // Kd is weighted if a vector has been passed by R. Else, a single numeric value has been passed.
  bool Weighted = (Weights.length() > 1);

  int d=0;
  for (int i=0; i < (x.length()-1); i++) {
    // Point j is a neighbor of i
    for (int j=i+1; j < x.length(); j++) {
      // i and j must be reference and neighbor
      if ((IsReferenceType[i] & IsNeighborType[j]) | (IsReferenceType[j] & IsNeighborType[i])) {
        // Calculate distance
        Dist[d] = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]));
        if (Weighted) {
          // if weighted, calculate the product of weights
        Weights[d] = PointWeight[i]*PointWeight[j];
        }
        d++;
      }
    }
  }
}

// [[Rcpp::export]]
void CountNbdKd(SEXP Rr, SEXP Rx, SEXP Ry, SEXP RWeight, SEXP RNbd, SEXP RIsReferenceType, SEXP RIsNeighborType) {
// Count the number of neighbors around each point for Kd, same as CountNbd but 
   // consider neighbors of interest only
   // do no attribute neighbors to each point, mix them
// Only ReferenceType points are considered. 
// The weights of NeighborType points are counted.

  //Distances
  NumericVector r(Rr);
  // x, y coordinates of points
  NumericVector x(Rx);
  NumericVector y(Ry);
  // Point weights
  NumericVector Weight(RWeight);
  // Matrix (single row) counting the number of neighbors. Modified by this routine.
  NumericMatrix Nbd(RNbd);
  // Boolean vectors describing reference and neighbor points
  IntegerVector IsReferenceType(RIsReferenceType);
  IntegerVector IsNeighborType(RIsNeighborType);

  double Distance2;
  double Nr = r.length();
  NumericVector r2 = r*r;
  int k = 0; 
  
  for (int i=0; i < (x.length()-1); i++) {
    // Consider reference type points
    if (IsReferenceType[i]) {
      // Point j is a neighbor of i. No neighbor is ignored.
      for (int j=i+1; j < x.length(); j++) {
        // Calculate squared distance
        Distance2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
        if (Distance2 <= r2[Nr-1]) {
          // Find the column of the matrix corresponding to the distance
          k = 0; 
          while (Distance2 > r2[k]) {
            k++;
          }
        } else {
          // Extra column for pairs far away
          k = Nr;
        }
        // The neighbor is a point of interest
        if (IsNeighborType[j]) {
          Nbd(0, k) += Weight[i]*Weight[j];
        }
        // j is a reference point
        if (IsReferenceType[j]) {
          // i is a point of interest around j
          if (IsNeighborType[i]) {
            Nbd(0, k) += Weight[i]*Weight[j];
          }
        }        
      }
    } else {
      // Point i is not a reference point
      for (int j=i+1; j < x.length(); j++) {
        // If point j is a reference point, it may be in its neighborhood
        if (IsReferenceType[j]) {
          // Calculate squared distance
          Distance2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
          if (Distance2 <= r2[Nr-1]) {
            // Find the column of the matrix corresponding to the distance
            k = 0; 
            while (Distance2 > r2[k]) {
              k++;
            }
          } else {
            // Extra column for pairs far away
            k = Nr;
          }
          // i is a point of interest around j
          if (IsNeighborType[i]) {
            Nbd(0, k) += Weight[i]*Weight[j];
          }
        }
      }
    }
  }
}
