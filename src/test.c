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
        //Builds an initial list of data points
        int n;
        for(n = 1;n<=SIZE;n++) {
            x[n] = n;
        }
        
        int Cos[SIZE];
        int Sin[SIZE];
        #ifdef USE_LIFT
        int RFactors[SIZE];
        #endif
        //1 2 3 4 5 6 7 8
        for (int j = 1; j <= (SIZE / 4);j++) {
            //1 5 9 13 17 21 25 29
            //This same angle generation was in original gg90()
            int angnum = 4 * ( j - 1 ) +1;
            int ang = angnum * PI / SIZE;
            //TODO: Cache/pre-calculate this
            Cos[j] = cos( ang );
            Sin[j] = sin( ang );
            #ifdef USE_LIFT
            //Reciprocal? Only in lifting. Dr. Mugler said it is "that constant"
            //I'm calling it the R Factor
            RFactors[j] = ( C[j] - 1 ) / S[j];
            #endif
            
            int c = Cos[j];
            int s = Sin[j];
            //0 2 4 6 8 10 12 14
            int cj = 2 * ( j - 1 );
            //31
            int m = SIZE - 1;

            
            //The tournament shuffle.
            //Takes first unused x[]
            //0 2 4 6 8 10 12 14
            temp[0] = x[cj];
            //Takes last unused x[]
            //31 29 27 25 23 21 19 17
            temp[1] = x[ ( m - cj ) ];
            //Takes first middle unused x[]
            //16 18 20 22 24 26 28 30
            temp[2] = x[ ( cj + ( SIZE / 2 ) ) ];
            //Second middle unused x[]
            //15 13 11 9 7 5 3 1
            temp[3] = x[ ( m - ( cj + (SIZE / 2) ) ) ];
            
            double v[4];
            #ifdef USE_LIFT
            lift(x[cj], x[ ( m - cj ) ], v[0], v[1], s, RFactors[j]);
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

//Rotates two points given in (x,y)
lift ( double x, double y, double xOut, double yOut, double sinValue, double RFactor ) {
    yOut = sinValue * (x - RFactor * y) - y;
    xOut = (x - RFactor * y) + RFactor * yOut;
}

lift ( double x[2], double out[2] double sinValue, double RFactor )
{
    out[1] = sinValue * (x[0] - RFactor * x[1]) - x[1];
    out[0] = (x[0] - RFactor * x[1]) + RFactor * out[1];
}

lift90sr ( double in[], double out[], double sinValues[SIZE], double RFactors[SIZE], int size )
{
    for (int i = 0; i < size; i+=2) {
        lift(in[i], in[i+1], &out[i], &out[i+1], sinValues[i], RFactors[i]);
    }
}

//Applies trig magic to every group of 4 in in[]
gg90 ( double in[], double out[] int size )
{
    for ( int i = 0; i < size; i+=4 ) {
        double angle = ( PI * i+1 ) / ( 2 * size );
        
        //Why aren't we caching this?
        double c = cos ( angle );
        double s = sin ( angle );
        
        out[i] = ( ( c * in[i] ) + ( s * in[i+1] ) );
        out[i+1] = ( s * in[i] ) - ( c * in[i+1] );
        out[i+2] = ( -s * in[i+2] ) + ( c * in[i+3] );
        out[i+3] = ( c * in[i+2] ) + ( s * in[i+3] );
        
    }
}