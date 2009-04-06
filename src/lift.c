#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#define PI 3.141592654
#define SQH 0.707106
#define size 32
#define b 4
double temp1[size];
double sum2[size];
double pout[size];
double poutb[size];
double loc1[size];
main ( int argc, char *argv[] )
{
    int L, m2, myrank, pcounter,hin, h5, h6, h7, h16, h8, J1, liftL, mrs, tag =
        1234, j4r[size],count, h1, h3, h4, numprocs, n = 0, i, j, cj, m, k = 1;
    double a1, a2, cf, sf, rf, ak, h2[10];
    MPI_Status status;
    double angnum, ang, mult, C[size], S[size], temp[size], temprs[size],
    v[size], S2[size], Rnew[size], Snew[size], bm, cm, p1[size], S4[size],
    R4[size], R2[size], temp2rs[size], p2[size], mytime, pair1[size],
    pair2[size], R[size], xtemp[size], Rold[size], temp2[size], Sold[size],
    temp3[size], temp4[size];
    double x[size], x1[size], c, s, r, s1, r1, p1a[size], p2a[size], p[size],
    xtemp1[size],fx1[size],fp[size];

    MPI_Init ( &argc, &argv );
    mytime = MPI_Wtime ();
    MPI_Comm_rank ( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size ( MPI_COMM_WORLD, &numprocs );

    if ( myrank == 0 ) { //procs 0 starts in here
        do {  //taking 8 values in x array
            x[n] = n + 1;
            n++;
        } while ( n < size );
        for ( j = 1; j <= ( n / 4 ); j++ ) {
            angnum = 4 * ( j - 1 ) + 1;
            mult = PI / ( n * 4 );
            ang = angnum * mult;
            C[j] = cos ( ang );
            S[j] = sin ( ang );
            R[j] = ( C[j] - 1 ) / S[j];
            //printf("s value is %lf,R %lf\n", S[j],R[j]);
        }
        for ( j = 1; j <= ( n / 4 ); j++ ) {
            cj = 2 * ( j - 1 );
            m = n - 1;
            //0,7,4,3 in the first iteration and 2 5 6 1 in the second iteration
            pair1[0] = x[cj];
            pair1[1] = x[ ( m - cj ) ];
            pair2[0] = x[ ( cj + ( n / 2 ) ) ]; //xtemp(3:4)
            pair2[1] = x[ ( m - ( cj + ( n / 2 ) ) ) ];

            lift ( pair1, S[j], R[j] );

            v[0] = pout[0];
            v[1] = pout[1];

            sumdiff ( pair2, b / 2 );

            pair2[0] = temp1[0] * SQH;
            pair2[1] = ( -1 * temp1[1] ) * SQH;

            lift ( pair2, S[j], R[j] );

            v[2] = pout[0];
            v[3] = pout[1];

            sumdiff ( v, b );

            p1[ ( 2 * j ) - 2] = temp1[0];
            p1[ ( 2 * j ) - 1] = temp1[1];
            p2[ ( 2 * j ) - 2] = temp1[2];
            p2[ ( 2 * j ) - 1] = temp1[3];
        }

        m = size / 2;
        J1 = 3;
        mrs = 4;

        for ( i = 0; i < mrs; i++ ) {
            temprs[i] = 1 + ( R[i + 1] * S[i + 1] );
            S2[i] = 2 * S[i + 1] * temprs[i];
            R2[i] = ( -1 * S[i + 1] ) / temprs[i];
            temp2rs[i] = 1 + ( R2[i] * S2[i] );
            S4[i] = 2 * S2[i] * temp2rs[i];
            R4[i] = ( -1 * S2[i] ) / temp2rs[i];
        }
        pcounter = 1;
        while ( pcounter <= 2 ) {
            if ( pcounter == 1 ) {
                for ( i = 0; i < n / 2; i++ ) {
                    p[i] = p1[i];
                }
            } else {
                m = size / 2;
                lift90sr ( p2, S4, R4, m );
                for ( i = 0; i < n / 2; i++ ) {
                    p[i] = poutb[i];
                }
            }
            liftL = 4;
            for ( i = 0; i < liftL; i++ ) {
                Rold[i] = R4[i];
                Sold[i] = S4[i];
            }
            for ( j = 1; j <= 2; j++ ) {
                m2 = pow ( 2, ( j - 1 ) );
                L = 16 / m2;
                for ( i = 0; i < liftL; i++ ) {
                    temprs[i] = 1 + Rold[i] * Sold[i];
                    Rnew[i] = -Sold[i] / temprs[i];
                    Snew[i] = 2 * Sold[i] * temprs[i];
                }
                for ( k = 1; k <= m2; k++ ) {
                    h1 = L * ( k - 1 );
                    h3 = ( L * k ) - 1;
                    h6 = L * ( k - 1 );
                    h16 = 0;

                    for ( h1; h1 <= h3; h1++ ) {
                        temp[h1] = p[h1];
                        h16++;
                    }
                    sumdiff2 ( temp, h6, h3 + 1 );
                    for ( i = 0; i < h16; i++ ) {
                        temp2[i] = sum2[i];
                    }
                    h4 = L / 2;
                    for ( i = 0, h4; h4 < L; h4++, i++ ) {
                        temp3[i] = temp2[h4]; //8-1a5
                    }

                    lift90sr ( temp3, Snew, Rnew, i + 1 );
                    for ( i = 0, h5 = L / 2; h5 < L; i++, h5++ ) {
                        temp2[h5] = poutb[i];
                    }

                    h7 = L * ( k - 1 );
                    h8 = ( L * k ) - 1;

                    for ( h7, i = 0; h7 <= h8; i++, h7++ ) {
                        p[h7] = temp2[i];
                    }
                }  //kloop ends here
                liftL = liftL / 2;
                for ( i = 0; i < 4; i++ ) {
                    Rold[i] = Rnew[i];
                    Sold[i] = Snew[i];
                }

            }   //end of for loop....
//   printf ("\n****************\n");
            //  for (i = 0; i < 16; i++)
            //  printf ("hi %lf\n", p[i]);

            j = J1;
            m2 = m / 4;
            L = 4;

            for ( k = 1; k <= m2; k++ ) {
                h1 = L * ( k - 1 );
                h3 = ( L * k ) - 1;
                h6 = L * ( k - 1 );
                h16 = 0;

                for ( h1; h1 <= h3; h1++ ) {
                    temp[h1] = p[h1];
                    h16++;
                    //      printf("temp is here%lf\n",temp[h1]);
                }
                sumdiff2 ( temp, h6, h3 + 1 );
                for ( i = 0; i < h16; i++ ) {
                    temp2[i] = sum2[i];
                    //printf("Temp2 is here %lf\n",temp2[i]);
                }
                h4 = L / 2;
                for ( i = 0, h4; h4 < L; h4++, i++ ) {
                    temp3[i] = temp2[h4];     //8-1a5
                }

                sumdiff2 ( temp3,0,i );
                for ( i = 0, h5 = L / 2; h5 < L; i++, h5++ ) {
                    temp2[h5] = sum2[i]*SQH;
                    //  printf("sumdiff result is %lf\n",temp2[h5]);
                }
                h7 = L * ( k - 1 );
                h8 = ( L * k ) - 1;

                for ( h7, i = 0; h7 <= h8; i++, h7++ ) {
                    p[h7] = temp2[i];
                }
//  printf("\n*****************\n");
                locations ( n,pcounter );


                j4r[k-1] = loc1[k];

                fx1[ j4r[k-1] ]         =  temp2[0];
                fx1[ ( n-1 )-j4r[k-1] ]   =  temp2[1];
                fx1[ ( n/2 )-1-j4r[k-1] ] =  temp2[2];
                fx1[  j4r[k-1]+ ( n/2 ) ]  =  temp2[3];
            }//k ends here
            pcounter++;  //have to end after the program
        }
        for ( i=0;i<32;i++ )
            printf ( "fx is %lf\n",fx1[i] );

        mytime = MPI_Wtime () - mytime;
        mytime = mytime * 1000000;
        printf ( "Timing from rank %d is %lfus.\n", myrank, mytime );
    } //end of rank 0
    MPI_Finalize ();
    return 0;
}
locations ( int n1,int pcounter )
{
    int dc =4;
    int mc = n1/16;
    int noterms = 1;
    int fcount = n1/8;
    int nparts = n1/8;
    double Js = log ( nparts ) /log ( 2 );
    int j1,nj,j2,nj2;

    loc1[1]= pcounter-1;

    for ( j1=1;j1<=Js;j1++ ) {
        for ( j2=1;j2<=noterms;j2++ ) {
            nj =  mc+1+ ( 2* ( j2-1 ) *mc );
            nj2=  loc1[1+2* ( j2-1 ) *mc];
            loc1[nj] = ( dc-1 )-nj2;
        }
        dc = 2*dc;
        mc = mc/2;
        noterms= 2*noterms;
        fcount = fcount/2;
    }
}

lift ( double x[size], double s, double R )
{
    double L1;
    L1 = x[0] - R * x[1];
    pout[1] = s * L1 - x[1];
    pout[0] = L1 + R * pout[1];
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
sumdiff2 ( double xy3[size], int r1, int r2 )
{
    int f1, r, i, j, h, k, l;
    double x1[size], x2[size], xy2[size];

    r = r2 - r1;
    f1 = ( r / 2 );

    for ( i = r1, k = 0; i < r2; i++, k++ )
        xy2[k] = xy3[i];

    for ( i = 0, k = 0; i < f1; i++, k++ )
        x1[k] = xy2[i];

    for ( h = 0, j = f1; j < r; j++, h++ )
        x2[h] = xy2[j];

    for ( k = 0; k < f1; k++ )
        sum2[k] = x1[k] + x2[k];

    for ( l = f1, k = 0; l < r; l++, k++ )
        sum2[l] = x1[k] - x2[k];
}
sumdiff ( double xy[size], int r )
{
    int f1, f2, i, j, h, k, l;
    double x1[size], x2[size];

    f1 = ( r / 2 );
    f2 = ( r / 2 ) + 1;

    for ( i = 0; i < f1; i++ )
        x1[i] = xy[i];

    for ( h = 0, j = f1; j < r; j++, h++ )
        x2[h] = xy[j];

    for ( k = 0; k < f1; k++ )
        temp1[k] = x1[k] + x2[k];

    for ( l = f1, k = 0; l < r; l++, k++ )
        temp1[l] = x1[k] - x2[k];
}
