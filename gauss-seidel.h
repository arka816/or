#include <bits/stdc++.h>
#define ITERATION_LIMIT 100
#define THRESHOLD 0.001

using namespace std;


bool check_convergence(vector<float> x, vector<float> last_x){
    // check if consecutive improvement is above threshold
    float max_change = 0;
    for(int i = 0; i < x.size(); i++){
        max_change = max(max_change, abs(x[i] - last_x[i]));
    }
    return max_change > THRESHOLD;
}

vector<float> gauss_seidel(vector<vector<float>> A, vector<float> b){
    int m = A.size(), n = A[0].size();

    // initial guess as all zero
    vector<float> x(n, 0);
    vector<float> last_x(n, 0);

    int i, j, sum, iter = 0;

    do{
        for(i = 0; i < n; i++){
            sum = b[i];
            for(j = 0; j < i; j++){
                sum -= A[i][j] * x[j];
            }
            for(j = i+1; j < n; j++){
                sum -= A[i][j] * last_x[j];
            }
            x[i] = sum / A[i][i];
        }
        // maintain a copy of the values from the last iteration
        for(i = 0; i < n; i++){
            last_x[i] = x[i];
        }
        iter ++;
    }while(check_convergence(x, last_x) && iter < ITERATION_LIMIT);

    return x;
}
