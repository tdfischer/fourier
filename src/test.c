/**
 * test.c
 * Copyright (c) 2009 Trever Fischer <wm161@wm161.net>
 * A merge-in-progress of lift.c and gg90.c
 */
#define SIZE 32

//Uncomment to use lift instead of gg90
//#define USE_LIFT

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int startTime = MPI_Wtime();
    
    int processorID;
    int worldSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &processorID);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    
    if (processorID == 0) {
        //Changed this from while(n<size) loop for readability
        int n;
        for(n = 1;n<=SIZE;n++) {
            x[n] = n;
        }
        
        int Cos[SIZE];
        int Sin[SIZE];
        #ifdef USE_LIFT
        int R[SIZE];
        #endif
        for (int j = 1; j <= (n / 4);j++) {
            int angnum = 4 * ( j - 1 ) +1;
            int mult = PI / ( n * 4 );
            int ang = angnum * mult;
            Cos[j] = cos( ang );
            Sin[j] = sin( ang );
            #ifdef USE_LIFT
            R[j] = ( C[j] - 1 ) / S[j]; //Reciprocal? Only in lifting. Dr. Mugler said it is "that constant"
            #endif
            
            int c = Cos[j];
            int s = Sin[j];
            int cj = 2 * ( j - 1 ); //Shows up in second loop in lifting
            int m = n - 1; //This too

            temp[0] = x[cj]; //And all this, but called "pair1/pair2" instead
            temp[1] = x[ ( m - cj ) ];
            temp[2] = x[ ( cj + ( n / 2 ) ) ];
            temp[3] = x[ ( m - ( cj + (n / 2) ) ) ];
            #ifdef USE_LIFT
            lift(temp, S[j], R[j]);
            #else
            
            #endif
        }
    }
}

void sumdiff ( double in[], double out[], int size )
{
    for (int i = 0; i < size / 2; i++ ) {
        //Split in[] into halves
        //Put first half plus second half into first half of return
        //And then the second half
        out[i] = in[i] + in[i*2];
        out[i*2] = in[i] - in[i*2];
    }
}

lift ( double x[2], double out[2] double sinValue, double RFactor )
{
    out[1] = sinValue * (x[0] - RFactor * x[1]) - x[1];
    out[0] = (x[0] - RFactor * x[1]) + RFactor * out[1];
}

lift90sr ( double pin[size], double Sm[size], double Rm[size], int m )
{
    int j, mj, i, N;
    double s[size], kemp[size], r[size], pout1[size], pout2[size], temp[size];
    
    for ( j = 1; j <= m / 4; j++ ) {
        s[j - 1] = Sm[j - 1];
        r[j - 1] = Rm[j - 1];
        mj = 2 * j - 1;
        kemp[0] = pin[mj - 1];
        kemp[1] = pin[mj];
        
        lift ( kemp, s[j - 1], r[j - 1] );
        pout1[mj - 1] = pout[0];
        pout1[mj] = pout[1];
        
        temp[0] = pin[ ( m / 2 ) + mj - 1];
        temp[1] = pin[ ( m / 2 ) + mj];
        
        lift ( temp, s[j - 1], r[j - 1] );
        pout2[mj - 1] = ( -1 ) * pout[1];
        pout2[mj] = pout[0];
    }
    for ( i = 0; i < m / 2; i++ )
        poutb[i] = pout1[i];
    
    for ( i = m / 2, j = 0; i < m; i++, j++ )
        poutb[i] = pout2[j];
}

gg90 ( double g2[size], int m2 )
{
    int j2 = 1, mj2, i, p; //m2 conatins teh length of array and g2 has the array elements of p2
    double Cm2[size], Sm2[size], c, s, angles2;
    
    
    for ( j2 = 1; j2 <= ( m2 / 4 ); j2++ ) {
        mj2 = 4 * ( j2 - 1 ) + 1;
        angles2 = ( PI * mj2 ) / ( 2 * m2 );
        Cm2[j2] = cos ( angles2 );
        Sm2[j2] = sin ( angles2 );
        
        c = Cm2[j2];
        s = Sm2[j2];
        
        p = 4 * ( j2 - 1 ) + 1;
        
        pout[mj2 - 1] = ( ( c * g2[p - 1] ) + ( s * g2[p] ) ); //overwriting the values check that
        pout[mj2] = ( s * g2[p - 1] ) - ( c * g2[p] );
        pout[mj2 + 1] = ( -s * g2[p + 1] ) + ( c * g2[p + 2] );
        pout[mj2 + 2] = ( c * g2[p + 1] ) + ( s * g2[p + 2] );
        
    }
}