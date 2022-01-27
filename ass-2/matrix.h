#include <bits/stdc++.h>
#include <assert.h>

using namespace std;

void print1d(vector<float> b){
    for(float v : b) cout<<v<<" ";
    cout<<endl;
}

void print2d(vector<vector<float>> A){
    for(vector<float> row : A){
        for(float v : row){
            cout<<v<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

vector<vector<float>> multiplyMatMat(vector<vector<float>> A, vector<vector<float>> B){
    int mA = A.size(), nA = A[0].size();
    int mB = B.size(), nB = B[0].size();

    assert(nA == mB);

    vector<vector<float>> prod(mA, vector<float> (nB, 0));
    int i, j, k;
    float sum;

    for(i = 0; i < mA; i++){
        for(j = 0; j < nB; j++){
            sum = 0.0;
            for(k = 0; k < nA; k++){
                sum += (float)(A[i][k] * B[k][j]);
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
    int i, j;
    float sum;

    for(i = 0; i < mA; i++){
        sum = 0;
        for(j = 0; j < nA; j++){
            sum += A[i][j] * b[j];
        }
        prod[i] = sum;
    }

    return prod;
}

vector<float> multiplyVecMat(vector<float> a, vector<vector<float>> B){
    int na = a.size();
    int mB = B.size(), nB = B[0].size();

    assert(na == mB);

    vector<float> prod(nB, 0);
    float sum;

    for(int i = 0; i < nB; i++){
        sum = 0;
        for(int j = 0; j < na; j++){
            sum += a[j] * B[j][i];
        }
        prod[i] = sum;
    }

    return prod;
}

vector<float> mulScalar(vector<float> a, float s){
    vector<float> prod(a);
    for(int i = 0; i < a.size(); i++) prod[i] = a[i] * s;
    return prod;
}

vector<float> addVec(vector<float> a, vector<float> b){
    int m = a.size();
    int n = b.size();

    assert(m == n);

    vector<float> sum(m, 0);

    for(int i = 0; i < m; i++){
        sum[i] = a[i] + b[i];
    }

    return sum;
}

float dot(vector<float> A, vector<float> B){
    float sum = 0;
    for(int i = 0; i < A.size(); i++){
        sum += A[i] * B[i];
    }
    return sum;
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

vector<vector<float>> transposeSq(vector<vector<float>> mat){
    int m = mat.size(), n = mat[0].size();
    assert(m == n);

    float temp;
    for(int i = 0; i < m; i++){
        for(int j = 0; j < i; j++){
            temp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = temp;
        }
    }

    return mat;
}

vector<vector<float>> inverse(vector<vector<float>> mat){
    int m = mat.size(), n = mat[0].size();
    assert(m == n);
    float det = determinant(mat);

    if(det == 0){
        return {};
    }

    vector<vector<float>> inv(m, vector<float> (n));

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            inv[i][j] = cofactor(mat, i, j) / det;
        }
    }

    return transposeSq(inv);
}
