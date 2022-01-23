/*
    lambda function - syntax used for capturing variables
    [&] : capture all external variable by reference 
    [=] : capture all external variable by value 
    [a, &b] : capture a by value and b by reference
*/

#include <bits/stdc++.h>
#include "gauss-seidel.h"
#include "matrix.h"

using namespace std;

vector<vector<int>> indices_list;

void combinationUtil(vector<int> curr, int start, int end, int index, int r){
    if (index == r){
        indices_list.push_back(curr);
        return;
    }

    for (int i = start; i <= end && end - i + 1 >= r - index; i++)
    {
        curr[index] = i;
        combinationUtil(curr, i+1, end, index+1, r);
    }
}

void print1d(vector<float> b){
    for(float v : b) cout<<v<<" ";
}

void print2d(vector<vector<float>> A){
    for(vector<float> row : A){
        for(float v : row){
            cout<<v<<" ";
        }
        cout<<endl;
    }
}

vector<vector<float>> bfs(vector<vector<float>> A, vector<float> b){
    // solves a lpp for basic feasible solutions written in standard form

    int m = A.size(), n = A[0].size();
    int i, j;
    vector<vector<float>> solns;

    // insert coefficients for slack variables
    for(i = 0; i < m; i++){
        A[i].resize(m+n, 0);
        A[i][n+i] = 1;
    }

    // generate combination vectors of non basic variable indices
    vector<int> curr(n, 0);
    combinationUtil(curr, 0, m+n-1, 0, n);

    for(vector<int> indices : indices_list){
        // set some variables as non-basic variables as zero
        // and solve the resultant system of equations
        vector<vector<float>> A_mod(A);

        reverse(indices.begin(), indices.end());
        // remove the corresponding columns from A
        for(int index : indices){
            // for(vector<float>& row : A){
            //     row.erase(next(row.begin(), index));
            // }
            for_each(A_mod.begin(), A_mod.end(), [&](vector<float>& row) -> void {
                row.erase(next(row.begin(), index));
            });
        }

        // cout<<"A modified: "<<endl;
        // print2d(A_mod);
        // cout<<endl;

        // solve 
        // vector<float> soln_mod = gauss_seidel(A_mod, b);
        vector<vector<float>> A_mod_inv = inverse(A_mod);

        // cout<<"A mod inverse: "<<endl;
        // print2d(A_mod_inv);
        // cout<<endl;

        if(A_mod_inv.size() == 0){
            cout<<"underdetermined system of equations"<<endl;
            continue;
        }
        vector<float> soln_mod = multiplyMatVec(A_mod_inv, b);

        
        reverse(indices.begin(), indices.end());
        // re-insert non-basic variables
        for(int index : indices){
            soln_mod.insert(soln_mod.begin() + index, 0);
        }

        // check if basic solution is feasible
        bool flag = true;
        for(int v : soln_mod){
            if(v < 0){
                flag = false;
                break;
            }
        }
        if(flag) solns.push_back(soln_mod);
    }

    return solns;
}

int main(){
    int n, m;
    vector<vector<float>> A;
    vector<float> b;
    vector<float> c;

    cout<<"Enter number of decision variables: "<<endl;;
    cin>>n;
    cout<<"Enter number of functional constraints: "<<endl;;
    cin>>m;
    cout<<"Enter the inequation matrix A (of dimensions "<<m<<" x "<<n<<" ): "<<endl;;
    float val;
    for(int i = 0; i < m; i++){
        vector<float> row(n, 0);
        A.push_back(row);
        for(int j = 0; j < n; j++){
            cin >> val;
            A[i][j] = val;
        }
    }
    cout<<"Enter the inequation vector b: "<<endl;;
    for(int i = 0; i < m; i++){
        cin >> val;
        b.push_back(val);
    }
    cout<<"Enter objective vector: "<<endl;
    for(int i = 0; i < m; i++){
        cin >> val;
        c.push_back(val);
    }

    vector<vector<float>> solns = bfs(A, b);

    cout<<endl<<"basic feasible solutions: "<<endl<<endl;
    for(vector<float> soln : solns){
        cout<<"solution: "<<"[";
        for(float v : soln){
            cout<<v<<", ";
        }
        cout<<"] "<<endl;
        cout<<"value of objective function: "<<dot(c, soln)<<endl<<endl;
    }

    return 0;
}