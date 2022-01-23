#include <bits/stdc++.h>
#include <algorithm>
#include "matrix.h"

using namespace std;

void print1d(vector<int> a){
    for(int v : a) cout<<v<<" ";
    cout<<endl;
}


float simplex(vector<vector<float>> A, vector<float> b, vector<float> c){
    int m = A.size(), n = A[0].size();

    int i, j;

    // create X_B : list of basic variables
    vector<int> bv;
    for(i = n; i < m+n; i++){
        bv.push_back(i);
    }

    // modify A
    // insert coefficients for slack variables
    for(i = 0; i < m; i++){
        A[i].resize(m+n, 0);
        A[i][n+i] = 1;
    }

    // modify c
    // insert coefficients for slack variables
    c.resize(m+n, 0);

    vector<vector<float>> B;
    vector<vector<float>> B_inv;
    vector<float> c_B;
    float opt = 0;

    vector<vector<float>> A_curr(A);
    vector<float> b_curr(b);
    vector<float> c_curr(c);

    int iterations = 0;

    while(true){
        // update basic variables
        // decide entering variable
        cout<<"iteration: "<<++iterations<<endl<<endl;

        int entering = 0;
        float min_c = c_curr[0];

        for(i = 1; i < m+n; i++){
            if(c_curr[i] < min_c){
                min_c = c_curr[i];
                entering = i;
            }
        }

        cout<<"entering variable: "<<entering<<endl;

        // decide leaving variable
        // min ratio test
        int leaving = -1;
        float min_ratio = FLT_MAX, ratio;

        for(j = 0; j < m; j++){
            if(b_curr[j] >= 0 && A_curr[j][leaving] > 0){
                ratio = (float)(b_curr[j] / A_curr[j][leaving]);
                if(ratio < min_ratio){
                    min_ratio = ratio;
                    leaving = j;
                }
            }
        }

        cout<<"leaving variable: "<<bv[leaving]<<endl;

        bv.erase(next(bv.begin(), leaving));
        bv.push_back(entering);

        cout<<"basic variables: ";
        print1d(bv);
        cout<<endl;

        B.clear();
        // construct B
        for(i = 0; i < m; i++){
            vector<float> row;
            for(int j : bv){
                row.push_back(A_curr[i][j]);
            }
            B.push_back(row);
        }

        c_B.clear();
        // construct c_B
        for(int j : bv){
            c_B.push_back(-1 * c_curr[j]);
        }

        B_inv = inverse(B);
        A_curr = multiplyMatMat(B_inv, A);
        b_curr = multiplyMatVec(B_inv, b);
        c_curr = addVec(c, multiplyVecMat(c_B, A_curr));
        opt = dot(c_B, b_curr);

        cout<<"A (curr): "<<endl;
        print2d(A_curr);
        cout<<endl;
        cout<<"b (curr): "<<endl;
        print1d(b_curr);
        cout<<endl;
        cout<<"c (curr): "<<endl;
        print1d(c_curr);
        cout<<endl;
        cout<<"optimum: "<<opt<<"\n\n\n\n";

        // if updated profit coefficients are all non-negative stop
        bool flag = true;
        for(i = 0; i < m+n; i++){
            if(c_curr[i] < 0){
                flag=false;
                break;
            }
        }
        if(flag) break;
    }

    return opt;
}

int main(){
    int n, m;
    vector<vector<float>> A{{3, 5}, {0, 1}, {8, 5}};
    vector<float> b{150, 20, 300};
    vector<float> c{-50, -40};

    // cout<<"Enter number of decision variables: "<<endl;;
    // cin>>n;
    // cout<<"Enter number of functional constraints: "<<endl;;
    // cin>>m;
    // cout<<"Enter the inequation matrix A (of dimensions "<<m<<" x "<<n<<" ): "<<endl;;
    // float val;
    // for(int i = 0; i < m; i++){
    //     vector<float> row(n, 0);
    //     A.push_back(row);
    //     for(int j = 0; j < n; j++){
    //         cin >> val;
    //         A[i][j] = val;
    //     }
    // }
    // cout<<"Enter the inequation vector b: "<<endl;;
    // for(int i = 0; i < m; i++){
    //     cin >> val;
    //     b.push_back(val);
    // }
    // cout<<"Enter objective vector: "<<endl;
    // for(int i = 0; i < m; i++){
    //     cin >> val;
    //     c.push_back(val);
    // }

    float optimum = simplex(A, b, c);
    cout<<optimum<<endl;
}