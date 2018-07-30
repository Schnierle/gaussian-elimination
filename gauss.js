var abs = Math.abs;

function array_fill(i, n, v) {
    var a = [];
    for (; i < n; i++) {
        a.push(v);
    }
    return a;
}

/**
 * Gaussian elimination
 * @param  array A matrix
 * @param  array x vector
 * @return array x solution vector
 */
function gauss(A, x) {

      //Index variables
      var i, k, j;

      var A_copy=[];
      var x_copy=[];

      //Copy A
      for (var i=0; i<A.length; ++i)
      {
        A_copy[i] = A[i].slice();
      }

      //Copy b
      for (var i=0; i<x.length; ++i)
      {
        x_copy[i] = x[i];
      }



    // Just make a single matrix
    for (i=0; i < A_copy.length; i++) {
        A_copy[i].push(x_copy[i]);
    }
    var n = A_copy.length;

    for (i=0; i < n; i++) {
        // Search for maximum in this column
        var maxEl = abs(A_copy[i][i]),
            maxRow = i;
        for (k=i+1; k < n; k++) {
            if (abs(A_copy[k][i]) > maxEl) {
                maxEl = abs(A_copy[k][i]);
                maxRow = k;
            }
        }


        // Swap maximum row with current row (column by column)
        for (k=i; k < n+1; k++) {
            var tmp = A_copy[maxRow][k];
            A_copy[maxRow][k] = A_copy[i][k];
            A_copy[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (k=i+1; k < n; k++) {
            var c = -A_copy[k][i]/A_copy[i][i];
            for (j=i; j < n+1; j++) {
                if (i===j) {
                    A_copy[k][j] = 0;
                } else {
                    A_copy[k][j] += c * A_copy[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    x_copy = array_fill(0, n, 0);
    for (i=n-1; i > -1; i--) {
        x_copy[i] = A_copy[i][n]/A_copy[i][i];
        for (k=i-1; k > -1; k--) {
            A_copy[k][n] -= A_copy[k][i] * x_copy[i];
        }
    }

    return x_copy;
}

module.exports = gauss;
