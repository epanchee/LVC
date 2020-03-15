/** @file
 *  @brief We search tips of a scroll wave here.
 *
 *  We do it using Kirsten ten Tusscher algorithm implemented by Hans Dierckx.
 *
/** @file
 *  @brief We search tips of a scroll wave here.
 *
 *  We do it using Kirsten ten Tusscher algorithm implemented by Hans Dierckx.
 *
 */

#ifndef TIPHANS_H_INCLUDED
#define TIPHANS_H_INCLUDED

/// we have nodes (0...nx, 0...ny)

/// maximal number of tips
#define maxntips 100

void track_tipline(int nx, int ny, double **var1, double **var2,
                   double* tipvals, double tipdata_[maxntips+1][2], int *tipsfound_);

#endif // TIPHANS_H_INCLUDED
