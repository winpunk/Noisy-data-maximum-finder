#include<stdio.h>
#include<math.h>

#define NZEROS 4
#define NPOLES 4
#define GAIN   2.400388646e+03



static float xv[NZEROS+1], yv[NPOLES+1];

void gaussEliminationLS(int m, int n, double a[m][n], double x[n - 1]) {
    int i, j, k;
    for (i = 0; i < m - 1; i++) {
        //Partial Pivoting
        for (k = i + 1; k < m; k++) {
            //If diagonal element(absolute vallue) is smaller than any of the terms below it
            if (fabs(a[i][i]) < fabs(a[k][i])) {
                //Swap the rows
                for (j = 0; j < n; j++) {
                    double temp;
                    temp = a[i][j];
                    a[i][j] = a[k][j];
                    a[k][j] = temp;
                }
            }
        }
        //Begin Gauss Elimination
        for (k = i + 1; k < m; k++) {
            double  term = a[k][i] / a[i][i];
            for (j = 0; j < n; j++) {
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }

    }
    //Begin Back-substitution
    for (i = m - 1; i >= 0; i--) {
        x[i] = a[i][n - 1];
        for (j = i + 1; j < n - 1; j++) {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }

}

double* Fit(int dataPoints, double x[dataPoints], double y[dataPoints])
{
    int N = dataPoints;//8;
    //degree of polynomial
    int n = 2;


    int i, j;



    // an array of size 2*n+1 for storing N, Sig xi, Sig xi^2, ...., etc. which are the independent components of the normal matrix
    double X[2 * n + 1];
    for (i = 0; i <= 2 * n; i++) {
        X[i] = 0;
        for (j = 0; j < N; j++) {
            X[i] = X[i] + pow(x[j], i);
        }
    }
    //the normal augmented matrix
    double B[n + 1][n + 2];
    // rhs
    double Y[n + 1];
    for (i = 0; i <= n; i++) {
        Y[i] = 0;
        for (j = 0; j < N; j++) {
            Y[i] = Y[i] + pow(x[j], i) * y[j];
        }
    }
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
            B[i][j] = X[i + j];
        }
    }
    for (i = 0; i <= n; i++) {
        B[i][n + 1] = Y[i];
    }
    double A[n + 1];
    printf("The polynomial fit is given by the equation:\n");
    //printMatrix(n + 1, n + 2, B);
    gaussEliminationLS(n + 1, n + 2, B, A);
    for (i = 0; i <= n; i++) {
        printf("%lfx^%d+", A[i], i);
    }
    return A;

}



static double* ButterworthArrayFilter(int size, double* indata, double deltaTimeinsec, double CutOff) {
    //if (indata == null) return null;
    if (CutOff == 0) return indata;

    double Samplingrate = 1 / deltaTimeinsec;
    long dF2 = size - 1;        // The data range is set with dF2
    double Dat2[dF2 + 4]; // Array with 4 extra points front and back
    double* data = indata; // Ptr., changes passed data

    // Copy indata to Dat2
    long r, s, t, p;

    for (r = 0; r < dF2; r++) {
        Dat2[2 + r] = indata[r];
    }
    Dat2[1] = Dat2[0] = indata[0];
    Dat2[dF2 + 3] = Dat2[dF2 + 2] = indata[dF2];

    const double pi = 3.14159265358979;
    double wc = tan(CutOff * pi / Samplingrate);
    double k1 = 1.414213562 * wc; // Sqrt(2) * wc
    double k2 = wc * wc;
    double a = k2 / (1 + k1 + k2);
    double b = 2 * a;
    double c = a;
    double k3 = b / k2;
    double d = -2 * a + k3;
    double e = 1 - (2 * a) - k3;

    // RECURSIVE TRIGGERS - ENABLE filter is performed (first, last points constant)
    double DatYt[dF2 + 4];
    DatYt[1] = DatYt[0] = indata[0];
    for (s = 2; s < dF2 + 2; s++) {
        DatYt[s] = a * Dat2[s] + b * Dat2[s - 1] + c * Dat2[s - 2]
                   + d * DatYt[s - 1] + e * DatYt[s - 2];
    }
    DatYt[dF2 + 3] = DatYt[dF2 + 2] = DatYt[dF2 + 1];

    // FORWARD filter
    double DatZt[dF2 + 2];
    DatZt[dF2] = DatYt[dF2 + 2];
    DatZt[dF2 + 1] = DatYt[dF2 + 3];
    for (t = -dF2 + 1; t <= 0; t++) {
        DatZt[-t] = a * DatYt[-t + 2] + b * DatYt[-t + 3] + c * DatYt[-t + 4]
                    + d * DatZt[-t + 1] + e * DatZt[-t + 2];
    }

    // Calculated points copied for return
    for (p = 0; p < dF2; p++) {
        data[p] = DatZt[p];
    }

    return data;
}

int FindMaxLocation (int size, double array[size])
{
    int c, location;
    double maximum = 0;

    for (c = 1; c < size; c++)
  {
    if (array[c] > maximum)
    {
       maximum  = array[c];
       location = c+1;
    }
  }
    return location;
}

main() {
    //no. of data-points
    double x_data[40];
    //double y_filtered[110];
    double* y_filtered;
    double y_peak_zone[40];

    int max_location;
    int i;

    for(i = 0; i<40; i++)
    {
        x_data[i]= i+1;

    }

    double y_data[110] =
    {
       0.596824089	,
0.482264232	,
0.447449748	,
0.562550252	,
0.477735768	,
0.503175911	,
0.339039052	,
0.405491503	,
0.332696703	,
0.180814427	,
0.18	,
0.170403548	,
0.162169262	,
0.195434697	,
0.2303301	,
0.116977778	,
0.015491503	,
0.095975952	,
0.158526204	,
0.003227271	,
0.04015369	,
0.059369152	,
0.0709262	,
0.004865966	,
-0.038782025	,
0.01	,
0.101217975	,
0.054865966	,
-0.0590738	,
0.109369152	,
-0.00984631	,
0.033227271	,
0.068526204	,
0.035975952	,
0.175491503	,
0.096977778	,
0.1203301	,
0.065434697	,
0.132169262	,
0.300403548	,
0.17	,
0.360814427	,
0.282696703	,
0.445491503	,
0.369039052	,
0.353175911	,
0.357735768	,
0.542550252	,
0.487449748	,
0.582264232	,
0.486824089	,
0.720960948	,
0.694508497	,
0.787303297	,
0.709185573	,
0.75	,
0.839596452	,
0.727830738	,
0.934565303	,
0.9196699	,
0.843022222	,
0.954508497	,
0.874024048	,
0.971473796	,
0.896772729	,
0.92984631	,
1.080630848	,
0.9190738	,
0.995134034	,
0.938782025	,
1.1	,
0.988782025	,
0.965134034	,
1.0490738	,
0.960630848	,
1.01984631	,
0.936772729	,
1.041473796	,
0.964024048	,
0.824508497	,
0.803022222	,
0.9596699	,
0.784565303	,
0.737830738	,
0.859596452	,
0.8	,
0.689185573	,
0.587303297	,
0.644508497	,
0.720960948	,
0.496824089	,
0.512264232	,
0.527449748	,
0.382550252	,
0.517735768	,
0.403175911	,
0.419039052	,
0.255491503	,
0.372696703	,
0.180814427	,
0.16	,
0.220403548	,
0.172169262	,
0.225434697	,
0.1503301	,
0.146977778	,
0.195491503	,
0.065975952	,
0.068526204	,
0.043227271
};

    y_filtered = ButterworthArrayFilter(110, y_data, 0.01, 5);

    for(i = 0; i < 110; i++)
    {        
        printf("%f\n", y_filtered[i]);
    }

    max_location = FindMaxLocation(110, y_filtered);
    printf("%d\n", max_location);

    printf("\n\nPeak zone:\n");
    for(i = (max_location - 20); i < (max_location+20); i++)
    {
        y_peak_zone[i-(max_location - 20)] = y_filtered[i];
        printf("%f\n", y_peak_zone[i-(max_location - 20)]);
    }

    double *fit = Fit(40, x_data, y_peak_zone);
    printf("\n%f, %f, %f", fit[0], fit[1], fit[2]);
    double a = fit[0] ;
    double b = fit[1] ;
    double c = fit[2] ;

    double temp[400];
    for(i=0; i < 400; i++)
    {

        temp[i] = (double)(i*i*c/100 + b*i/10 + a);
        printf("\n%lf", temp[i]);
    }

    printf("\n%d",FindMaxLocation(400, temp));
}
