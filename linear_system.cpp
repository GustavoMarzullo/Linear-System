#include <iostream>
#include <vector>
#include <math.h>

using std::cin; using std::cout; using std::endl;
using std::vector;

typedef vector<vector<double>> matrix;
typedef vector<double> doublevec;

double detTriang(matrix A){
    int n=A.size();
    int p = 1; 
    int k;
    for(k = 0; k < n - 1; k++){
        double max = abs(A[k][k]);
        int maxIndex = k;
        int i;
        for(i = k + 1; i < n; i++){
            if(max < abs(A[i][k])){
                max = abs(A[i][k]);
                maxIndex = i;
            }
        }
        if(maxIndex != k){
            double temp;
            p = p * (-1);
            int j;
            for(j = 0; j < n; j++){
                temp = A[k][j];
                A[k][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }
        }
        if(A[k][k] == 0){
            return 0;
        }else{
            int m;
            for(m = k + 1; m < n; m++){
                double F = (-1) * A[m][k]/A[k][k];
                A[m][k] = 0;
                int l;
                for(l = k + 1; l < n; l ++){
                    A[m][l] = A[m][l] + F * A[k][l];
                }
            }
        }
    }
    double det = 1;
    int q;
    for(q = 0; q < n; q++){
        det = det * A[q][q];
    }
    return p * det;
}

matrix GetVariableMatrix(int n){
    int i,j;
    vector<vector<double>> A(n, vector<double> (n+1));
    printf("Enter the A matrix\n");
    for(int i=0;i<n;i++){
     for(int j=0;j<n;j++){
        scanf("%lf",&A[i][j]);
        }
    }
    printf("\n");
    return A;
}

doublevec GetConstantVector(int n){
    int i;
    doublevec b(n,0);
    printf("Enter the b matrix\n");
    for(i=0;i!=n;++i){
        scanf("%lf",&b[i]);
    }
    printf("\n");
    return b;
}

double _xi(matrix A, doublevec b, double det_A, int j){
    int n=A.size(),i;
    for(i=0;i!=n;++i){
        A[i][j]=b[i];
    }
    return detTriang(A)/det_A;
}

void cramer(matrix A, doublevec b, int n){
    doublevec results(n);
    double det_A=detTriang(A);
    for(int i=0;i!=n;++i){
        results[i]=_xi(A,b,det_A,i);
    }
    printf("Results:\n");
    for(int i=0;i!=n;++i){
        printf("x%d = %lf\n",i+1,results[i]);
    }
}

void print(const matrix &A){
    int i,j;
    int m,n;
    m=A.size();
    n=A[0].size();
    for(i=0;i!=m;++i){
        for(j=0;j!=n;++j){
            cout<<A[i][j]<<" ";
        }
        cout << endl; //new line 
    }
    printf("\n");
}

int main(){
    int n;
    printf("Linear system (A*x=b) solver\n");
    printf("How many variables does your system have?\n");
    scanf("%d",&n);
    printf("\n");
    matrix A;
    A=GetVariableMatrix(n);
    doublevec b;
    b=GetConstantVector(n);
    cramer(A,b,n);
    printf("\n");
    return 0;
}


