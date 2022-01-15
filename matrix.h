#include <bits/stdc++.h>
#include <assert.h>

using namespace std;

vector<vector<float>> multiplyMatMat(vector<vector<float>> A, vector<vector<float>> B){
    int mA = A.size(), nA = A[0].size();
    int mB = B.size(), nB = A[0].size();

    assert(nA == mB);

    vector<vector<float>> prod(mA, vector<float> (nB));
    int i, j, k, sum;

    for(i = 0; i < mA; i++){
        for(j = 0; j < nB; j++){
            sum = 0;
            for(k = 0; k < nA; k++){
                sum += A[i][k] + B[k][j];
            }
            prod[i][j] = sum;
        }
    }

    return prod;
}

vector<float> multiplyMatVec(vector<vector<float>> A, vector<float> b){
    int mA = A.size(), nA = A[0].size();
    int mB = b.size();

    assert(nA == mB);

    vector<float> prod(mA, 0);
    int i, j, sum;

    for(i = 0; i < mA; i++){
        sum = 0;
        for(j = 0; j < nA; j++){
            sum += A[i][j] * b[j];
        }
        prod[i] = sum;
    }

    return prod;
}

vector<vector<float>> minor(vector<vector<float>> mat, int col){
    vector<vector<float>> minor(mat);
    // erase first row
    minor.erase(minor.begin());
    // erase column
    for_each(minor.begin(), minor.end(), 
        [&](vector<float>& row) -> void {
            row.erase(next(row.begin(), col));
        }
    );
    return minor;
}

float determinant(vector<vector<float>> mat){
    int m = mat.size(), n = mat[0].size();
    assert(m == n);
    float det = 0;

    if(m == 1){
        return mat[0][0];
    }

    for(int i = 0; i < m; i++){
        det += (i % 2 == 0 ? 1 : -1) * determinant(minor(mat, i)) * mat[0][i];
    }
    return det;
}

float cofactor(vector<vector<float>> mat, int row, int col){
    vector<vector<float>> cofactor(mat);
    // erase row
    cofactor.erase(next(cofactor.begin(), row));
    // erase col
    for_each(cofactor.begin(), cofactor.end(), 
        [&](vector<float>& row) -> void {
            row.erase(next(row.begin(), col));
        }
    );
    return ((row + col) % 2 == 0 ? 1 : -1) * determinant(cofactor);
}

vector<vector<float>> inverse(vector<vector<float>> mat){
    int m = mat.size(), n = mat[0].size();
    assert(m == n);
    float det = determinant(mat);
    assert(det != 0);

    vector<vector<float>> inv(m, vector<float> (n));

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            inv[i][j] = cofactor(mat, i, j) / det;
        }
    }

    return inv;
}
